module Proof

open System.Runtime.InteropServices

/// sharing polynomial indexes with C++
type whichPoly = V = 0 | W = 1 | Y = 2
let V = whichPoly.V
let W = whichPoly.W
let Y = whichPoly.Y
let poly = function 0 -> V | 1 -> W | 2 -> Y | _ -> failwith "Bad poly" in

/// bank type, from the verifier's viewpoint
type bankType = RECOMPUTED = 0 | VERIFIED = 1 | BOTH = 2
let RECOMPUTED = bankType.RECOMPUTED
let VERIFIED   = bankType.VERIFIED
let BOTH       = bankType.BOTH

type whichKey = EK = 1 | VK = 2 | BOTH = 3


module KeyGen = 
  [<DllImport(@"FFLib.dll")>] extern void private key_delete (void*)
  [<DllImport(@"FFLib.dll")>] extern void private save_key (void*, void*, whichKey, string, bool)
  [<DllImport(@"FFLib.dll")>] extern void* private load_key (void*, whichKey, string)
  type k = // generated keys
    { value: nativeint } with  
    override k.Finalize() = () // key_delete(k.value) // disable to prevent multiple deallocations of reloaded keys
    static member load which filename = { value = load_key(Field.encoding, which, filename) }
    static member save which filename key simple_format = save_key(Field.encoding, key.value, which, filename, simple_format) 

  [<DllImport(@"FFLib.dll")>] extern void* private key_gen_init (void*, int, int, int[], bankType[])
  [<DllImport(@"FFLib.dll")>] extern void  private key_gen_add_poly_elt(void*, int, int, whichPoly, int, void*)
  [<DllImport(@"FFLib.dll")>] extern void  private key_gen_add_poly_elt_by_val(void*, int, int, whichPoly, int, int)
  [<DllImport(@"FFLib.dll")>] extern void* private key_gen_finalize (void* )
  type t =  
    { value: nativeint } // internal state, to be used linearly (KeyGenContext* in cpp); explicitly finalized below.
    
  // both arrays are indexed by static banks, indicating their size & bankType     
  let init (roots:int) (banks: int array, bts: bankType array) = 
    assert (banks.Length = bts.Length)
    { value = key_gen_init (Field.encoding, roots, banks.Length, banks, bts) }

  let inline add (state:t) b k p r (x:Field.t) = 
    key_gen_add_poly_elt (state.value,b,k,p,r,x.raw)

  let inline addI state b k p r v = 
    if Lib.vv then printfn "  bank%3d: %A[r%3d] = %A x_%d " b  p r v k
    if v <= bigint System.Int32.MaxValue && v >= bigint System.Int32.MinValue // && false to disable 
    then key_gen_add_poly_elt_by_val (state.value,b,k,p,r,int v)
    else add state b k p r (Field.t.ofBig v)

  let finalize (state:t): k =
    { value = key_gen_finalize(state.value) }


[<DllImport(@"FFLib.dll")>] extern void private proof_delete(void *)
[<DllImport(@"FFLib.dll")>] extern void private save_proofs(void *, int, void* [], string, bool)
[<DllImport(@"FFLib.dll")>] extern void* private load_proof(void *, string)

/// QAP proof: a few field elements (excluding bank commitments), wrapping newProof*
/// each proof file may contains a series of related proofs, but we use only single ones for now
type proof = 
    { value : nativeint } 
    override x.Finalize() = proof_delete x.value
    static member load (file:string) =  { value =  load_proof(Field.encoding, Lib.path + "proof-" + file + ".proof") }
    member x.Save mode file =
      let fullname = "proof-" + file
      match mode with
      | Lib.Simple -> save_proofs(Field.encoding, 1, [| x.value |], Lib.path + fullname + ".data", true)
      | Lib.Rich   -> save_proofs(Field.encoding, 1, [| x.value |], Lib.path + fullname + ".proof", false)

[<DllImport(@"FFLib.dll")>] extern void* private commit_init (void*, bool, int)


// Consider using the n-ary "commits" variant
[<DllImport(@"FFLib.dll")>] extern void private commit_print (void* )
[<DllImport(@"FFLib.dll")>] extern void* private load_commit (void*, void*, string)
[<DllImport(@"FFLib.dll")>] extern void private save_commit (void*, void*, string, bool)

[<DllImport(@"FFLib.dll")>] extern void private commit_nonzero_poly  (void*, whichPoly, int, void*, void* )
[<DllImport(@"FFLib.dll")>] extern void private commit_enc_poly      (void*, int, void* )
[<DllImport(@"FFLib.dll")>] extern void private commit_enc_poly_multi(void*, int[], int, int, void* )
[<DllImport(@"FFLib.dll")>] extern void private commit_finalize      (void*)
[<DllImport(@"FFLib.dll")>] extern void private commit_streamline    (void*)
[<DllImport(@"FFLib.dll")>] extern void private commit_delete        (void*)

/// Commitments
/// prover: also with partial polynomial sum on the other side
/// verifier: commitment to all values for some given key & bank index
type commitment = 
    { value: nativeint } with 
    override x.Finalize() =  commit_delete(x.value)
    member x.Save mode file =
      let fullname = "commit-" + file
      if Lib.v then printfn "saving %s" fullname 
      match mode with 
      | Lib.Simple -> save_commit(Field.encoding, x.value, Lib.path + fullname + ".data", true)
      | Lib.Rich   -> save_commit(Field.encoding, x.value, Lib.path + fullname + ".commit", false)
      
    static member load (k:KeyGen.k) file = 
      { value = load_commit(Field.encoding, k.value, Lib.path + "commit-" + file + ".commit") }

let zeros: commitment = { value = System.IntPtr.Zero }


let npoly = ref 0 // number of prover polynomial coefficients
let nenc  = ref 0 // number of prover exponentials

module Prove = 
  type key = EK of KeyGen.k 
  let loadKey file : key = EK(KeyGen.k.load whichKey.EK (Lib.path + file + ".ek"))
  let saveKey mode file (EK key) = 
    match mode with 
    | Lib.Simple -> KeyGen.k.save whichKey.EK (Lib.path + file + ".data") key true
    | Lib.Rich   -> KeyGen.k.save whichKey.EK (Lib.path + file + ".ek")   key false

  /// building commitments, one variable assignment at a time; then we finalize; then we save.
  let commitment (EK key) bankId : commitment = { value = commit_init(key.value, true, bankId) }

  // set a non-zero value at root --- doesn't matter which polynomial this is coming from.
  // we always accumulate dense polynomials (towards the global division)
  // we independently accumulate bank commitments, unless the variable is public  
  let inline add_poly (c:commitment) (p:whichPoly) r (x:Field.t) (v:Field.t) = 
    if Lib.vv then printfn " |  addpoly %A[%A] = %A * %A" p r v x  
    incr npoly;
    commit_nonzero_poly(c.value, p, r, x.raw, v.raw)
  let inline add (c:commitment) k (x:Field.t) = 
    if Lib.vv then printfn " |  addc    x%d   = %A" k x  
    incr nenc
    commit_enc_poly(c.value, k, x.raw)
  // batched, same as "for i in st..en-1 do add c k.[i] x 
  let inline add_multi (c:commitment) (ks: int[]) st en (x:Field.t) = 
    nenc := !nenc + en - st
    if Lib.vv then 
      printfn "add_multi %d..%d" st (en - 1)
      for i = st to en - 1 do
        printfn " |  addc    x%03d   = %A" ks.[i] x  
    commit_enc_poly_multi(c.value, ks, st, en, x.raw)
  
  // complete the commitment
  let finalize  (c:commitment) = commit_finalize(c.value)

  // make the commitment usable only for verification (discarding the polynomials)
  let streamline (c:commitment) = commit_streamline(c.value)
  
  [<DllImport(@"FFLib.dll")>] extern void* private prove (void*, void* []) 
  
  /// compute a proof for a collection of commitments (one for each bank)
  let proof (EK key) (cs: commitment array) : proof = { value = prove(key.value, Array.map (fun (c:commitment) -> c.value) cs) }

let venc  = ref 0 // number of verifier exponentials

module Verify =
  type key = VK of KeyGen.k 
  let loadKey file : key    = VK(KeyGen.k.load whichKey.VK (Lib.path + file + ".vk"))
  let saveKey mode file (VK key) = 
    match mode with 
    | Lib.Simple -> KeyGen.k.save whichKey.VK (Lib.path + file + ".data" ) key true
    | Lib.Rich   -> KeyGen.k.save whichKey.VK (Lib.path + file + ".vk")    key false

  /// building commitments, one variable assignment at a time; then finalize.
  let load (VK key) name = commitment.load key name 
  let recompute (VK key) bankId : commitment = { value = commit_init(key.value, false, bankId) }
  let inline add (c:commitment) k (x:Field.t) = 
    if Lib.vv then printfn " |  addc    x%d   = %A" k x  
    incr venc
    commit_enc_poly(c.value, k, x.raw)
  let inline add_multi(c:commitment) k st en (x:Field.t) = 
    venc := !venc + en - st
    commit_enc_poly_multi(c.value,k,st,en,x.raw)
  let finalize (c:commitment) : unit = commit_finalize(c.value)

  [<DllImport(@"FFLib.dll")>] extern int private verify_commit (void*, int, void*, bool)
  let committed (VK key) bankId (c:commitment) = 
    verify_commit(key.value, bankId, c.value, Lib.expect_to_fail) <> 0
  
  [<DllImport(@"FFLib.dll")>] extern int private verify_proof (void*, void*, int, void* [], bool)
  let proof (VK key) (pi: proof) (cs: commitment array) = 
    // printfn "cs = %A" cs 
    // Array.iter (fun c -> if c = zeros then printfn "zeros" else commit_print c) cs 
    let nbanks = cs.Length
    let mutable ok = true
    let cs' = Array.map (fun (c:commitment) -> c.value) cs
    for i = 1 to 1 do // just for measuring the average 
      ok <- ok && verify_proof(key.value, pi.value, nbanks, cs', Lib.expect_to_fail) <> 0
    ok





module Test =

  let basic() = 
      Lib.expect_to_fail <- true
      // equation x0 * z = x1 + 2t
      let r = 0
      // KEYGEN
      let s = KeyGen.init 1 ([| 2; 1; 1 |],[| RECOMPUTED; VERIFIED; BOTH |]) // 1 root, 2 IO variables, 1 local variable, 1 shared variable
      KeyGen.add s 0 0 V r Field.t.one //V: x0
      KeyGen.add s 0 1 Y r Field.t.one //Y: x1 
      KeyGen.add s 1 0 W r Field.t.one //W: z
      KeyGen.add s 2 0 Y r Field.t.two //Y: 2t
      let ek, vk = let k = KeyGen.finalize s in Prove.EK k, Verify.VK k
      // EVAL
      let x0 = Field.t.ofInt32 7
      let x1 = Field.t.ofInt32 13
      let z  = Field.t.ofInt32 3
      let t  = Field.t.ofInt32 4
      // PROVER
      let c0 = Prove.commitment ek 0  
      Prove.add_poly c0 V r Field.t.one x0 // duplicating  V: x0 so we must keep the IO polynomials for proving.
      Prove.add_poly c0 Y r Field.t.one x1 // duplicating  Y: x1
      let c1 = Prove.commitment ek 1
      Prove.add c1 0 z
      Prove.add_poly c1 W r Field.t.one z
      let c2 = Prove.commitment ek 2
      Prove.add c2 0 t
      Prove.add_poly c2 Y r Field.t.two t
      let pi = Prove.proof ek [| c0; c1; c2 |] 
      Prove.finalize c0
      Prove.finalize c1
      Prove.finalize c2
      c1.Save Lib.Rich "test"
      pi.Save Lib.Rich "test"
      let proof = (pi,c1,c2) // the communicated proof 
      // VERIFIER
      let c0 = Verify.recompute vk 0 
      Verify.add c0 0 x0
      Verify.add c0 1 x1
      Verify.finalize c0
      // the c0 commitment is locally computed, not cryptographically verifiable.
      printfn "Verify.committed vk 1                : %b" (Verify.committed vk 1 c1)
      printfn "Verify.committed vk 2                : %b" (Verify.committed vk 2 c2)
      printfn "Verify.committed (bad c, should fail): %b" (Verify.committed vk 2 c1)
      printfn "Verify.committed (bad c, should fail): %b" (Verify.committed vk 1 c2)
      printfn "Verify.proof                         : %b" (Verify.proof vk pi [| c0; c1; c2 |])
      begin 
        let c0 = Verify.recompute vk 0 
        Verify.add c0 0 x0
        Verify.add c0 1 (x1 + Field.t.one)
        Verify.finalize c0
        printfn "Verify.proof (bad input, should fail): %b" (Verify.proof vk pi [| c0; c1; c2 |])
      end
      printfn "Verify.proof (duplicate, should fail): %b" (Verify.proof vk pi [| c0; c0; c2 |])
      printfn "Verify.proof (swapped; ok?)          : %b" (Verify.proof vk pi [| c1; c2; c0 |]) 
      printfn "Verify.proof (missing, should fail)  : %b" (Verify.proof vk pi [| c0; zeros; zeros |])
      // unclear how to simply test bad verify.committed 
  
      //printfn "result: %b" (Verify.result vk [| c0; c1 |] pi)
      //printfn "result: %b" (Verify.result vk [| zero; c1 |] pi)
      //printfn "result: %b" (Verify.result vk [| c0; zero |] pi)
  
   
      //TODO: define keypair as a private type, with finalizer too 

  open QAP

  let mini_fg() =

    // further reduced from the fg test
    let banks = 3 
    let roots = 2
    let a = KeyGen.init roots ([| 1; 2; 2 |], [| RECOMPUTED; RECOMPUTED; VERIFIED |])
    let add p b k r (v:int) = KeyGen.addI a b k p r (bigint v) 

    // equations for f:  root   0:  (x1)(x1) = (x4)
    //                   root   1:  (1 + x2 + x4)(1) = (x3)

    add V 0 0  1  1
    add W 0 0  1  1
    add V 1 0  0  1
    add W 1 0  0  1
    add V 1 1  1  1
    add Y 2 0  1  1
    add V 2 1  1  1
    add Y 2 1  0  1

    let ek, vk = let k = KeyGen.finalize a in Prove.EK k, Verify.VK k
    
    // prove from local run [[|1|]; [|2; 3|]; [|8|]; [|4|]] 
    let cs = Array.create banks zeros 
    let add p b k r v x =
      Prove.add_poly cs.[b] p r (Field.t.ofInt32 v) (Field.t.ofInt32 x)  
    let add2 p b k r v x = // only for true banks; only once per (b,k)
      add p b k r v x 
      Prove.add cs.[b] k (Field.t.ofInt32 x)

    cs.[0] <- Prove.commitment ek 0 
    add  V 0 0  1  1  1
    add  W 0 0  1  1  1
    Prove.finalize cs.[0]

    cs.[1] <- Prove.commitment ek 1 
    add  V 1 0  0  1  2
    add  W 1 0  0  1  2
    add  V 1 1  1  1  3
    Prove.finalize cs.[1]

    cs.[2] <- Prove.commitment ek 2 
    add2 Y 2 0  1  1  8
    add2 V 2 1  1  1  4
    add  Y 2 1  0  1  4
    Prove.finalize cs.[2]

    let pi = Prove.proof ek cs

    // verify from the same local run
    let csv = Array.create banks zeros 
    let add p b k r v x = Verify.add csv.[b] k (Field.t.ofInt32 x) 

    csv.[0] <- Verify.recompute vk 0 
    add V 0 0  1  1  1
    Verify.finalize csv.[0]

    csv.[1] <- Verify.recompute vk 1 
    add V 1 0  0  1  2
    add V 1 1  1  1  3
    Verify.finalize csv.[1]

    if Verify.committed vk 2 cs.[2] then csv.[2] <- cs.[2] 

    printfn "Verify.proof %b" (Verify.proof vk pi csv) 

  // the "fg" mQAP
  let fg() =

    // from running Bank.keygen, then manually removed linking
    let banks = 4 (*7*)
    let s = 3 (*6*)
    let roots = 3
    let link = true // test both!
    let a = KeyGen.init roots ([| 1; 2; 2; (* 1; 1; 2; *) 1  |], [| RECOMPUTED; RECOMPUTED; VERIFIED; (* RECOMPUTED; RECOMPUTED; VERIFIED; *) VERIFIED |])
    let add p b k r (v:int) = KeyGen.addI a b k p r (bigint v) 

    add V 0 0  1  1
    add W 0 0  1  1
    if link then add W 0 0  2  1

    add V 1 0  0  1
    add W 1 0  0  1
    add V 1 1  1  1
    
    if link then add V 2 0  2  1
    add Y 2 0  1  1
    add V 2 1  1  1
    add Y 2 1  0  1

    (*
    add V 3 0  1 -1
    add W 3 0  1  1
    if link then add W 3 0  2  1

    add Y 4 0  1  1

    add V 5 0  0  1
    add W 5 0  0  1
    if link then add V 5 0  2  1
    add V 5 1  1  1
    add Y 5 1  0  1
    *)

    add Y s 0  2  1

    let ek, vk = let k = KeyGen.finalize a in Prove.EK k, Verify.VK k

    // prove from local run [[|1|]; [|2; 3|]; [|8|]; [|4|]] 
    let cs = Array.create banks zeros 
    let add p b k r v x =
      Prove.add_poly cs.[b] p r (Field.t.ofInt32 v) (Field.t.ofInt32 x)  
    let add2 p b k r v x = // only for true banks; only once per (b,k)
      add p b k r v x 
      Prove.add cs.[b] k (Field.t.ofInt32 x)

    cs.[0] <- Prove.commitment ek 0 
    add  V 0 0  1  1  1
    add  W 0 0  1  1  1
    if link then add  W 0 0  2  1  1
    Prove.finalize cs.[0]

    cs.[1] <- Prove.commitment ek 1 
    add  V 1 0  0  1  2
    add  W 1 0  0  1  2
    add  V 1 1  1  1  3
    Prove.finalize cs.[1]

    cs.[2] <- Prove.commitment ek 2 
    if link then add  V 2 0  2  1  8
    add2 Y 2 0  1  1  8    
    add2 V 2 1  1  1  4
    add  Y 2 1  0  1  4
    Prove.finalize cs.[2]

    if link then 
      cs.[s] <- Prove.commitment ek s
      add2 Y s 0  2  1  8
      Prove.finalize cs.[s]

    let pi = Prove.proof ek cs

    // verify from the same local run
    let csv = Array.create banks zeros 
    let add p b k r v x = Verify.add csv.[b] k (Field.t.ofInt32 x) 

    csv.[0] <- Verify.recompute vk 0 
    add V 0 0  1  1  1
    Verify.finalize csv.[0]

    csv.[1] <- Verify.recompute vk 1 
    add V 1 0  0  1  2
    add V 1 1  1  1  3
    Verify.finalize csv.[1]

    if         Verify.committed vk 2 cs.[2] then csv.[2] <- cs.[2] 
    if link && Verify.committed vk s cs.[s] then csv.[s] <- cs.[s] 

    printfn "Verify.proof %b" (Verify.proof vk pi csv) 
 
  // what we test on each sample below:
  let fflib roots inputs outputs inters =
    let unit  = [| Field.t.one |]
    let io    = Array.map (Field.t.ofInt32) (Array.append (Array.ofList inputs) (Array.ofList outputs))
    let local = Array.map Field.t.ofInt32 (Array.ofList inters)

    printfn "%A %A %A" unit io local
    let a = KeyGen.init roots ([| 1; Array.length io; Array.length local; 2 |],[| RECOMPUTED; RECOMPUTED; VERIFIED; RECOMPUTED |]) 

    let k = ref 0
    let add b n =
       let w = getWire !k
       incr k
       Map.iter (KeyGen.addI a b n V) w.V  
       Map.iter (KeyGen.addI a b n W) w.W
       Map.iter (KeyGen.addI a b n Y) w.Y  
    add 0 0 
    Array.iteri (fun i _ -> add 1 i) io
    Array.iteri (fun i _ -> add 2 i) local 
    k := 0; add 3 0; add 3 1 // doesn't matter what those are
    let ek, vk = let k = KeyGen.finalize a in Prove.EK k, Verify.VK k

    k := 0
    let add' both c n v =
       let w = getWire !k
       incr k
       Map.iter (fun r x -> Prove.add_poly c V r (Field.t.ofBig x) v) w.V
       Map.iter (fun r x -> Prove.add_poly c W r (Field.t.ofBig x) v) w.W
       Map.iter (fun r x -> Prove.add_poly c Y r (Field.t.ofBig x) v) w.Y
       if both then Prove.add c n v
    let cs = Array.append (Array.init 3 (Prove.commitment ek)) [| zeros |] 
    Array.iteri (add' false cs.[0]) unit
    Array.iteri (add' false cs.[1]) io
    Array.iteri (add' true  cs.[2]) local
    Prove.finalize cs.[0]
    Prove.finalize cs.[1]
    Prove.finalize cs.[2]
    let pi = Prove.proof ek cs

    cs.[0] <- Verify.recompute vk 0
    cs.[1] <- Verify.recompute vk 1
    Array.iteri (fun i v -> Verify.add cs.[0] i v) unit
    Array.iteri (fun i v -> Verify.add cs.[1] i v) io
    Verify.finalize cs.[0]
    Verify.finalize cs.[1]
    printfn "Verify commitment: %A" (Verify.committed vk 2 cs.[2])
    printfn "Verify proof     : %A" (Verify.proof vk pi cs)
    
(* older:
  // what we test on each sample below:
  let fflib_old roots inputs outputs inters =
    let inputs  = Field.ints2field inputs
    let outputs = Field.ints2field outputs
    let inters  = Field.ints2field inters
    let ostart = inputs.size + inters.size + 1
    let qap = Native.pack()
    print()
    let keys = generateKeys roots QAP.wire.Count qap inputs.size outputs.size
    let proof = prove keys.pub qap keys.targetRoots inputs outputs inters 
    if not (verify keys proof qap inputs outputs ostart false)            then failwith "Verification should succeed!"
    let inputs' = Field.ints2field (List.init inputs.size (fun i -> i % 2))
    let outputs' = Field.ints2field (List.init outputs.size (fun i -> i % 2))
    if verify keys proof qap inputs outputs' ostart true    then failwith "Verification on bad output should fail!"
    if verify keys proof qap inputs' outputs ostart true then failwith "Verification on bad input should fail!"
*)

  let add() = 
    wire <-
        let r0 = newRoot()
        let r1 = newRoot()
        printfn "r0 = %d, r1 = %d" r0 r1
        let p = wire.[0] 
        let single r (v:int) = Map.empty.Add(r, bigint v)
        [| 
           {p with                   W = single 0 9;                  } 
           {p with V = single 0 1;                                    } 
           {p with                   W = single 0 1;                  } 
           {p with V = single 1 1;                                    } 
           {p with                   W = single 1 1;                  } 
           {p with                                     Y = single 0 1 } 
           {p with                                     Y = single 1 1 } |]
    let c1,c2,c3,c4 = 1,2,3,4
    let c5 = c1 * (c2+9) // equation 0  
    let c6 = c3 * c4     // equation 1
    let roots = 2
    fflib roots [ c1;c2;c3;c4 ] [ c5;c6 ] []

  let layered() = 
    wire <-
        let r0 = newRoot()
        let r1 = newRoot()
        let p = wire.[0] 
        let single r (v:int) = Map.empty.Add(r, bigint v)
        [|  p
            {p with V = single 0 1;                                       } 
            {p with                        W = single 0 1;                } 
            {p with                        W = single 1 1;                } 
            {p with V = single 1 1;                        Y = single 0 1 } 
            {p with                                        Y = single 1 1 } |]
            
    let c1,c2,c3 = 1,2,3
    let c4 = c1 * c2   // equation 0  
    let c5 = c4 * c3   // equation 1
    let inputs  = [ c1;c2;c3 ] // without c0!
    let outputs = [ c5 ]
    let inters  = [ c4 ]    
    let roots = 2
    let ostart = 5
    fflib roots inputs outputs inters

  let simple() =
    init() 
    wire <- 
      let r0 = newRoot()
      let r1 = newRoot()
      let p = wire.[0] 
      let single r (v:int) = Map.empty.Add(r, bigint v)
      [|  {p with V = single 0 1                                 } 
          {p with V = single 0 1;                                } 
          {p with                 W = single 0 1;                } 
          {p with V = single 1 1;                                } 
          {p with                 W = single 1 1;                } 
          {p with                                 Y = single 0 1 } 
          {p with                                 Y = single 1 1 } |]
    let c1,c2,c3,c4 = 1,2,3,4
    let c5 = (1 + c1) * c2   // equation 0  
    let c6 = c3 * c4   // equation 1
    let inputs  = [ c1;c2;c3;c4 ] 
    let outputs = [ c5 ]
    let inters  = [ c6 ]    
    let roots = 2
    fflib roots inputs outputs inters

  let broken() = // aligned to trivial.c
    init()
    wire <- 
        let r0 = newRoot()
        let r1 = newRoot()
        let p = wire.[0] 
        let single r (v:int) = Map.empty.Add(r, bigint v)
        [| {p with V = single r1 11; W = single r1 1;                }
           {p with V = single r0 7;                                  } 
           {p with                  W = single r0 1;                 } 
           {p with V = single r1 -1                                  } 
           {p with V = single r1 1;                  Y = single r0 1 } |]
    let c1,c2 = 3,5
    let c3 = 7 * c1 * c2 // equation 0  
    let c4 = 11 + c3     // coded as equation 1: (11 + c3 - c4)*1 = 0
    let inputs  = [ c1;c2 ] 
    let outputs = [ c4 ]
    let inters  = [ c3 ]  
    let roots = 2
    fflib roots inputs outputs inters 

  let poly() = // as above, except the coding of the last equation    
    init()
    wire <- 
        let r0 = newRoot()
        let r1 = newRoot()
        let p = wire.[0] 
        let single r (v:int) = Map.empty.Add(r, bigint v)
        [|  {p with V = single r1 11; W = single r1 1;                 }
            {p with V = single r0 7 ;                                  } 
            {p with                   W = single r0 1;                 } 
            {p with                                    Y = single r1 1 } 
            {p with V = single r1 1;                   Y = single r0 1 } |]
    let c1,c2 = 3,5
    let c3 = 7 * c1 * c2  // equation 0  
    let c4 = 11 + c3      // equation 1
    let inputs  = [ c1;c2 ] 
    let outputs = [ c4 ]
    let inters  = [ c3 ]  
    let roots = 2
    fflib roots inputs outputs inters 

  let poly_old() = 
    init()    
    wire <- 
        let r0 = newRoot()
        let r1 = newRoot()
        let p = wire.[0] 
        let single r (v:int) = Map.empty.Add(r, bigint v)
        [| {p with V = single r1 3; W = single r1 1;                 }
           {p with V = single r0 1;                                  } 
           {p with                  W = single r0 1;                 } 
           {p with                                   Y = single r1 1 } 
           {p with V = single r1 2;                  Y = single r0 1 } |]
    let c1,c2 = 3,5
    let c4 = c1 * c2        // equation 0  
    let c3 = (3 + 2*c4) * 1 // equation 1
    let inputs  = [ c1;c2 ] // without c0!
    let outputs = [ c3 ]
    let inters  = [ c4 ]  
    let roots = 2
    fflib roots inputs outputs inters 
