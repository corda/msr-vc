
module Pinoccioq

module QapReader = // for interop with compcert back-end

  let next_uint (h: System.IO.FileStream) : uint32 =
      let b3 = uint32(h.ReadByte()) in
      let b2 = uint32(h.ReadByte()) in
      let b1 = uint32(h.ReadByte()) in
      let b0 = uint32(h.ReadByte()) in
      let r = (b3 <<< 0o30) ||| (b2 <<< 0o20) ||| (b1 <<< 0o10) ||| b0
      r
  let next_int (h: System.IO.FileStream) = int (next_uint h)

  let next_bigint (h: System.IO.FileStream) =
      let mutable r = QAP.zero in
      for i = 0 to 7 do
        r <- (r <<< 32) + bigint (next_uint h) 
      r

  let read
    (* we read QAPS with interleaved banks; so each bank is an array of wires, in no particular order *)
    (* for now:  bank 0 is public "c0" 
                 bank 1 is for all IOs, a single public bank
                 bank 2 is for the locals *)
    (filename: string)
    (cb_init: int -> int array array -> 'state)               (* number of roots, wires per bank *)
    (cb_eqn: int -> int -> int -> bigint -> 'state -> 'state) (* polynomial, root, wire, coefficient *)
    : 'state
   =
    let h = System.IO.File.OpenRead(filename) in
    let magic = next_int h in
    if (magic <> 0x51415030) then failwithf "Bad magic number: %x" magic
    let nroots = next_int h;
    let nbanks = next_int h;
    let bank_wires = Array.create nbanks [||] in
    for b = 0 to nbanks - 1 do
      let sz = next_int h in
      bank_wires.[b] <- Array.create sz 0;
      for w = 0 to sz - 1 do
        bank_wires.[b].[w] <- next_int h - 1 //variable, reindexed from 0 instead of 1
      done
    done;
    let state = ref (cb_init nroots bank_wires) in
    for p = 0 to 2 do
      let s = next_int h in
      for i = 0 to s - 1 do
        let t = next_int h in
        for j = 0 to t - 1 do
          let r = next_int h - 1 in
          let w = next_int h - 1 in
          let k = next_bigint h in
          state := cb_eqn p r w k !state
        done
      done
    done;
    h.Close()
    !state

  // Renumber the wires to group them by bank.
  // This function computes a renumbering table indexed by old values.
  // Values in the table are pairs (wire, bank).
  let renumber (wires: int array array) : (int * int) array=
    let size = Array.fold (fun n a -> n + Array.length a) 0 wires in
    let fwd = Array.create size (0, 0) in
    Array.fold (fun i a ->
      Array.fold (fun n w ->
        fwd.[w] <- (n, i)
        n + 1
      ) 0 a |> (ignore : int -> unit);
      i + 1
    ) 0 wires |> (ignore : int -> unit);
    fwd

  type t = {
    nroots : int
  ; bank_sizes : int array
  ; bank_types : Proof.bankType array
  ; renumber : (int * int) array
  ; data : Map<int, bigint> array array array 
  }
   
  type v = Field.t array array

  let doit (filename: string) : t =
    read filename (fun nr bk ->
      let fwd = renumber bk in
      let bs = Array.map Array.length bk in
      { nroots = nr
      ; renumber = fwd
      ; bank_sizes = bs
      ; bank_types = [| Proof.RECOMPUTED; Proof.VERIFIED|]
      ; data = Array.init 2 (fun b -> Array.init bs.[b] (fun w -> Array.create 3 Map.empty)) }
    ) (fun p r w k a ->
      let (w, b) = a.renumber.[w] in
      a.data.[b].[w].[p] <- Map.add r k a.data.[b].[w].[p];
      a
    )

  let solution (qap:t) (filename: string) = 
    let v = Array.map2 (fun typ size -> Array.create size Field.t.zero) qap.bank_types qap.bank_sizes
    let h = System.IO.File.OpenRead(filename) in
    //let magic = next_int h in
    //if (magic <> 0x51415030) then failwithf "Bad magic number: %x" magic
    Array.iter (fun (i,b) -> v.[b].[i] <- Field.t.ofBig (next_bigint h)) qap.renumber
    h.Close()
    v


let run file datafile =
  printfn "Loading QAP from %s" file
  let qap = QapReader.doit file in
  let witnesses = QapReader.solution qap datafile
  if Lib.v then printfn "QAP = %A\n witnesses = %A\n" qap witnesses
  let keygen() = 
    let kg = Proof.KeyGen.init qap.nroots (qap.bank_sizes, qap.bank_types) in
    Array.iteri (fun b ->
      Array.iteri (fun k ->
        Array.iteri (fun p ->
          Map.iter (fun r v ->
            Proof.KeyGen.addI kg b k (Proof.poly p) r v)))) qap.data
    Proof.KeyGen.finalize kg

  let N = 1
  let rec repeat n f x = if n = 1 then f x else ignore(f x); repeat(n-1) f x  
  let k = repeat N (Lib.time "key generation" keygen) () 

  // Worker commitments
  let ek = Proof.Prove.EK k
  let cs = 
    Array.mapi (fun b ws ->
      let ws = witnesses.[b]
      let commit() = 
        let c = Proof.Prove.commitment ek b
        if qap.bank_types.[b] = Proof.VERIFIED then 
          Array.iteri (fun k witness -> Proof.Prove.add c k witness) ws
        Array.iteri (fun k prvs ->
          let x = ws.[k]
          Array.iteri (fun p ->
            Map.iter (fun r v ->  
              Proof.Prove.add_poly c (Proof.poly p) r x (Field.t.ofBig v)
          )) prvs )
          qap.data.[b]
        Proof.Prove.finalize(c)
        c 
      repeat N (Lib.time (sprintf "commitment %d generation" b) commit) ())
      witnesses

  // Worker proof generation
  let pi = repeat N (Lib.time "proof generation" (Proof.Prove.proof ek)) cs

  let vk = Proof.Verify.VK k
  // Verifier commitments 
  let csv = 
    Array.mapi2 (fun b ws c -> 
      if qap.bank_types.[b] = Proof.VERIFIED then 
        let ok = repeat N (Lib.time (sprintf "commitment %d verification" b) (Proof.Verify.committed vk b)) c
        assert(ok);
        c
      else
        let recompute() = 
          let c = Proof.Verify.recompute vk b
          Array.iteri (fun k witness -> Proof.Verify.add c k witness) ws
          Proof.Verify.finalize c
          c
        repeat N (Lib.time (sprintf "commitment %d recomputation" b) recompute) ())
      witnesses cs

  // Verifier proof checking
  let valid = repeat N (Lib.time "proof verification" (Proof.Verify.proof vk pi)) csv
  assert(valid)
