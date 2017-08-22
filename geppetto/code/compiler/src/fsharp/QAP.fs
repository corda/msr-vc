module QAP

// This module provides functions for constructing QAPs.
// Their representation as sparse polynomials is internal.

/// variable (wire) indexes
type k = int 

/// root indexes
type r = int 

/// interpreter representation for all integers at compile-time
type integer = bigint

// it is worth sharing common integer constants
let zero   = bigint.Zero
let one    = bigint 1
let two    = bigint 2
let minus1 = -one

/// integer intervals min..max such that -p < min <= max < p && max - min < p 
module Range =
  [<StructuredFormatDisplay("{AsString}")>]
  type r = { min: integer; max: integer } 
  with 
    member x.AsString = sprintf "[%A..%A]" x.min x.max

  let inline point (x:integer) = { min = x; max = x }

  let inline uint b: r = { min = zero; max = (one <<< b) - one } 
  let inline int b : r = { min = -one <<< (b-1); max = (one <<< (b-1)) - one } 
  let uint8  = uint  8
  let uint16 = uint 16
  let uint32 = uint 32

  let minField = one - Field.p 
  let maxField = Field.p - one
  let full = { min = zero; max = maxField }
  let inline fit (r:r) = minField <= r.min && r.max <= maxField 

  let bit = { min = zero; max = one }
  let zero = point zero
  let one = point one

  let inline neg x = { min = - x.max; max = - x.min }   
  let inline shift x n = { min = x.min + n; max = x.max + n }
  let inline inter x y = { min = max x.min y.min; max = min x.max y.max } 
  let inline within x r = r.min <= x && x <= r.max

  let inline (+) x y = { min = x.min + y.min; max = x.max + y.max }   
  let inline (-) x y = x + neg y
  let inline mulc x (c:integer)  = 
    if c.Sign = 1 
    then { min = x.min * c; max = x.max * c} 
    else { min = x.max * c; max = x.min * c} 
  let inline (*) x y = 
    let a = Array.sort [| x.min * y.min; x.min * y.max; x.max * y.min; x.max * y.max |]
    { min = a.[0]; max = a.[3] } 

    
module Semantics = 
  type r = Range.r

  /// C integers and field element semantics to (eventually) enforce
  /// note that the actual range may exceed it, e.g. as we delay truncations
  [<StructuredFormatDisplay("{AsString}")>]
  type t = 
  | Elt of r option // field semantics, possibly with precise range
  | Int of int * r  // binary semantics, of width n, with precise range (possibly larger than 2^n) 
  with member x.AsString = match x with 
                             | Elt None     -> "elt"
                             | Elt (Some r) -> sprintf "elt%A" r
                             | Int(n,r)     -> sprintf "int%d%A" n r

  let int1_zero = Int(1,Range.zero)
  let int1_one  = Int(1,Range.one)
  let int1 =  Int( 1, Range.bit)
  let int8   = Int( 8, Range.int  8)
  let int16  = Int(16, Range.int 16)
  let int32  = Int(32, Range.int 32) 
  let int64  = Int(64, Range.int 64) 

  let uint8  = Int( 8, Range.uint8) 
  let uint16 = Int(16, Range.uint16) 
  let uint32 = Int(32, Range.uint32) 

  let elt = Elt None
  let elt0 = Elt (Some Range.zero)
  let elt1 = Elt (Some Range.one)
  let eltbit = Elt (Some Range.bit)

  let range x = 
    match x with 
    | Elt(Some r) | Int(_,r) -> r
    | Elt None               -> Range.full

  let inline lax r = if Range.fit r then Elt(Some r) else elt
  let inline lax2 (f: r -> r -> r) o o' = 
    match o, o' with
    | None, _ | _, None -> elt
    | Some r0, Some r1  -> lax(f r0 r1) 
    
  let inline strict n r = 
    if Range.fit r then Int(n,r) 
    else  failwith "field overflow"
  let inline strict2 n n' r = 
    if Range.fit r then 
      if n = n' then Int(n,r) 
      else 
        if Lib.v then printfn "warning: implicit extension from int%d to int%d" (min n n') (max n n')
        Int(max n n',r) 
    else  failwithf "field overflow on int%d: %A" n r

  let neg x = 
    match x with 
    | Int(n,r)    -> Int(n, Range.neg r)
    | Elt(Some r) -> Elt(Some (Range.neg r))
    | _           -> elt
  let (+) x y = 
    match x, y with  
      | Int(n,r), Int(n',r') -> strict2 n n' (Range.(+) r r') 
      | Elt o, Elt o'        -> lax2 Range.(+) o o'   
      | Int(_,r), Elt o
      | Elt o, Int(_,r)      ->
         // we tolerate implicit coercions int --> elem, e.g. for conditionals
         if Lib.v then printfn "warning: mixed semantics %A + %A" x y 
         lax2 Range.(+) o (Some r)

  let (*) x y = 
    match x, y with  
      | Int(n,r), Int(n',r')             -> strict2 n n' (Range.(*) r r') 
      | _, _  when x = elt0 || y = elt0  -> elt0 // tolerating int 0 --> elt coercionse  
      | Int(_,r), _ when r = Range.bit   -> y    // tolerating int 0..1 --> elt coercions
      | Elt o, Elt o'                    -> lax2 Range.(*) o o'  
      | Int(_,r), Elt o
      | Elt o, Int(_,r)      ->
         // we tolerate implicit coercions int --> elem, e.g. for conditionals
         if Lib.v then printfn "warning: mixed semantics %A * %A" x y 
         lax2 Range.(*) o (Some r)

  let (/) x y = 
    match x, y with 
      | Elt _, Elt _ -> elt 
      | _            -> failwithf "mixed semantics %A + %A" x y
  let muln x (c:integer) =
    match x with 
      | Elt(Some r) -> if c.IsZero then elt0 else lax (Range.mulc r c)
      | Int(n,r)    -> if c.IsZero then Int(n,Range.zero) else strict n (Range.mulc r c) 
      | _           -> if c.IsZero then elt0 else elt

  let shift x n = 
    match x with 
    | Int(b,r)    -> strict b (Range.shift r n)
    | Elt(Some r) -> lax (Range.shift r n)
    | _           -> elt
  let inter x y = 
    match x, y with  
    | Int(n,r), Int(n',r') when n = n' -> strict n (Range.inter r r') 
    | Elt o, Elt o'                    -> lax2 Range.inter o o'   
    | _                                -> failwithf "mixed semantics %A + %A" x y
  let bound x r' = 
    match x with  
    | Int(n,r)    -> strict n (Range.inter r r') 
    | Elt(Some r) -> Elt(Some (Range.inter r r'))
    | Elt None    -> lax r'

  // turning bits into their proper semantics
  // we expect r to be 0..0, 0..1, or 1..1
  // we could check the range of s is wider
  let coerce_bit s s1 = 
    match s, s1 with 
    | Elt _   , Int(1,r) -> if   r = Range.bit  then eltbit 
                            elif r = Range.zero then elt0 
                            elif r = Range.one  then elt1
                            else failwith "Elt(Some r)"
    | Int(n,_), Int(1,r) -> Int(n,r) 
    | _                  -> failwithf "not a bit: %A" s1


// ------------------------------------------- splitting integers into bit vectors

let log2 (n: bigint) = // inefficient; will do for now.
    let mutable m = n
    let mutable l = 0
    assert(not m.IsZero)
    while not m.IsOne do
      m <- m >>> 1
      l <- l + 1
    l      

type sign = Signed | Unsigned
type plan = int * sign * integer // b,s,h
/// plan a succinct b-bit decomposition for a given range:
/// either unsigned:   h * 2^b <= ... <  (h + 1)   * 2^b 
/// or signed: (h - 1/2) * 2^b <= ... <  (h + 1/2) * 2^b 

let planSplit (r: Range.r) : plan = 
    let (b,s,h) as result =     
      let b = 
        if r.max = r.min 
        then 0 
        else log2 (r.max - r.min) + 1    // minimal number of bits needed 
      let m = r.min >>> (b - 1)          // m * 2^(b-1) <= r.min
      if r.max < ((m + two) <<< (b - 1)) // r.max < m * 2^(b-1) + 2^b ? 
      then
        if m.IsEven then b, Unsigned, m / two 
                    else b, Signed, (m + one) / two 
      else
        let b = b + 1
        let m = r.min >>> (b - 1)  
        if m.IsEven then b, Unsigned, m / two 
                    else b, Signed, (m + one) / two 

    if b > 254 then // should be 253
      printfn "planSplit: too many bits in %A" r
    if Lib.v then 
      if s = Unsigned 
      then printfn "Unsigned(%d,%A) : %A  <= %A <= v <= %A < %A" b h (h <<< b) r.min r.max ((h + one) <<< b)   
      else printfn "Signed(%d,%A) : %A  <= %A <= v <= %A < %A" b h ((two*h - one) <<< b - 1) r.min r.max ((two*h + one) <<< b - 1)
    result

let rec split n (v:integer) =
  let r = ref v 
  Array.init 
    n
    (fun i -> let x = if (!r).IsEven then 0 else 1 in r := !r >>> 1; x)

let merge (bs: int array) = // assuming unsigned semnatics
  let r = ref zero
  for i = bs.Length-1 downto 0 do
    r := two * !r + if bs.[i] = 1 then one else zero
  !r

/// splitting algorithm for integers (also used as a spec for linear)
/// returns both the b bits actually split and the resulting n bits
let splitInteger semantics (b,s,h) (v:integer) = 
    if Lib.v && not (Range.within v (Semantics.range semantics)) then 
      printfn "warning: %A not in %A" v semantics
    match s with
    | Unsigned  -> 
        let v' = v - (h <<< b)
        let bs = split b v'
        // printfn "split %A: %10A %3d <= %3d" semantics s n b
        let bits = 
          match semantics with
          | Semantics.Elt _    -> bs // keep all bits
          | Semantics.Int(n,_) -> if n <= b then Array.sub bs 0 n // immediately truncate
                                  else let hs = split (n - b) h in Array.append bs hs 
        bs, bits

    | Signed ->
        // We have a signed n-bit integer represented as a element mod p.
        let v' = v - (h <<< b) + (one <<< (b - 1)) 
        let v'' = v' % Field.p 
        let bs = split b v''
        let c = bs.[b-1] // 1 iff 2^(b-2) <= v' 
        match semantics with
        | Semantics.Elt _ -> failwithf "unexpected elt" 
        | Semantics.Int(n,_) -> 
          let bits = 
            if n <= b 
            then 
              //printfn "splitInteger %10A %3d <= %3d" s n b
              Array.sub bs 0 n 
            else
              let hs0 = split (n - b) (h - one) 
              let hs1 = split (n - b) h
              let hs = Array.map2 (fun h0 h1 -> h0 + c * (h1 - h0)) hs0 hs1
              Array.append bs hs
          if b <= n then bits.[b-1] <- 1 - c // flipping the sign bit  
          bs, bits

 // ------------------------------------------- hint log
        
module Log = 
  /// runtime hints, passed from generation to evaluation, logically part of the evaluation key.
  /// the record and hint must be synchronized.
  /// the order is significant; beware of nested calls logging other hints
  /// GEN:  init ; record* ; save 
  /// EVAL: load ; hint*
  type hint = 
    | Split of Semantics.t * bool * plan // What's the semantics? Shall I split? What is the representation?  
    | Mlt1 | Mlt0          // Shall I record the result of this multiplication? Cheaper than Mlt of bool
    | Zerop of bool        // Shall I zerop (using 2 variables)?
    | Truncate of bool     // Shall I truncate?
    | Output of bool       // Shall I assign? 
    | Nroot of int         // How many equations so far?   
    | IsBit of bool        // Is this int always in 0..1?    
    | Empty
  type log = Log of hint array

  let mutable private recorded = [||] 
  let mutable private loaded   = [||] 
  let private pos = ref 0
  let size = 1 <<< 20 // initial number and increment for recorded hints, ~ number of roots

  let init() = 
    loaded   <- [||]
    recorded <- Array.create size Empty
    pos := 0;

  let verbose h = () //printfn "record[%7d] : %A" !pos h
    
  let record h = 
    verbose h
    if !pos = recorded.Length then recorded <- Lib.resize recorded (recorded.Length + size) 
    recorded.[!pos] <- h; incr pos
  let save() =  
    let r = Log(Array.sub recorded 0 !pos)
    recorded <- [||]
    if Lib.v then printf "#hints = %A\n" !pos
    r

  let clear() = 
    loaded <- [||]
    recorded <- [||]

  let load (Log a) = pos := 0; loaded <- a
  let hint() = 
    let h = loaded.[!pos]
    verbose h 
    incr pos; h  

// ------------------------------------------- building QAPs with sparse polynomials

/// We record every non-zero value p(root[r])  
type poly = Map<r,integer> 
let p0 : poly = Map.empty
let p1 r = p0.Add(r,one)

type fpoly = (int * string)[]
let pformat (p:poly) : fpoly = Array.map (fun (r,i:integer) -> r, i.ToString()) (Map.toArray p)
let pparse  (f:fpoly):  poly = Map.ofArray (Array.map (fun (r,s:string) -> r, bigint.Parse s) f)

/// Partitioning wires; C++ relies instead on their numbering
/// (not used much for now)
type bank =
  | One               // special IO always equal to 1
  | Shared of string  // named shared bank, including I/O
  | Local             // local wires for intermediate values; can be freely eliminated

type wire = { mutable V: poly
              mutable W: poly
              mutable Y: poly 
              bank: bank 
              }
// we could use p indexed by V = 0 | W = 1 | Y = 2
// but then beware of arrays becoming mutable. 
 
let no_wire: wire[] = [||] 

// We keep some global mutable state (read-only outside this module)
let mutable nWire : k = 0
let mutable nRoot : r = 0
let mutable wire: wire array = no_wire

let free k = 
  let w = wire.[k]
  w.V = Map.empty && w.W = Map.empty && w.Y = Map.empty

// a (single) QAP interface consists of a series of banks,
// each one an array of a base type (for now, native | field) * a range [a QAP invariant]
// how to represent partiality? Do we care?
// do we need to truncate from one QAP to another? 

// a (single) QAP consists of 
// - a number of roots (its degrees)
// - a number of local wires (each with a range)
// - three sparse polynomials p[wire] 

let inline setV r k v = let x = wire.[k] in x.V <- x.V.Add(r,v)  // allocating 60% of managed memory :(
let inline setW r k v = let x = wire.[k] in x.W <- x.W.Add(r,v)
let inline setY r k v = let x = wire.[k] in x.Y <- x.Y.Add(r,v)

let newRoot _ = let r = nRoot in nRoot <- r + 1; r
let newWire w = let k = nWire in  
                if wire.Length = nWire then wire <- Lib.resize wire (nWire + nWire/3 + (1 <<< 16))
                nWire <- k + 1 
                wire.[k] <- w 
                k

open System.Collections
let splitCache : Hashtable = Hashtable()

let init() =
  // also clear any cache! 
  nWire <- 0;
  nRoot <- 0;
  wire  <- no_wire
  Log.init()
  splitCache.Clear()
  // the first wire (0); we'll keep extending its polynomials with the equation constants
  let k0 = newWire { V = p0; W = p0; Y = p0; bank = One }
  assert(k0 = 0)

let clear() = 
  nWire <- 0;
  nRoot <- 0;
  wire  <- no_wire
  Log.clear()
  splitCache.Clear()

type t = { nWire: k; nRoot: r; wire: wire[]; log: Log.log } 

let dummy = { nWire = 0; nRoot = 0; wire = no_wire; log = Log.Log [||] }

let nCoeff (qap:t)  = 
  let mutable v = 0
  let mutable w = 0
  let mutable y = 0

  for k = 0 to qap.nWire - 1 do 
    let x = qap.wire.[k] 
    v <- v + x.V.Count
    w <- w + x.W.Count
    y <- y + x.Y.Count
  (v,w,y)

let save() = 
  let q =    { nWire = nWire
               nRoot = nRoot
               wire  = Lib.resize wire nWire 
               log   = Log.save() }
  clear(); q

let load s = 
  nWire <- s.nWire
  nRoot <- s.nRoot
  if Lib.check then wire  <- s.wire
  Log.load s.log

let wformat w = w.bank, pformat w.V, pformat w.W, pformat w.Y
let wparse (b,v,w,y) = { bank = b; V = pparse v; W = pparse w; Y = pparse y }

type fwire = (bank * fpoly * fpoly * fpoly)[]
type f = k * r * fwire * Log.log
let format (qap:t) : f = qap.nWire, qap.nRoot, Array.map wformat qap.wire, qap.log 
let parse (k,r,w,l) = { nWire = k; nRoot = r; wire = Array.map wparse w; log = l } 

let getWire k = wire.[k] 

// these two are internal to QAP
let wire_iter f = for i = 0 to nWire - 1 do f i wire.[i]
let wire_fold f a =
  let mutable a = a
  for i = 0 to nWire - 1 do
    a <- f a i wire.[i]
  a

// this one is used on saved QAPs, e.g. in Bank
let wire_foldi f a (ws: wire[]) =
  let mutable a = a
  for i = 0 to ws.Length - 1 do
    a <- f a i ws.[i]
  a

/// Symbolic values are sparse linear combinations of wires 
/// with a target integer semantics and an approximate range
/// (invariant: integer coefficients are non-null) 
type linear = Map<k,integer> 
type value  = Semantics.t * linear 
let valueSemantics ((s,_):value) = s 

let linear1 k = Map.empty.Add(k, one) 

let get (a:linear) k             = match a.TryFind k with Some n -> n | None -> zero
let set (a:linear) k (v:integer) = if v.IsZero then a.Remove k else a.Add(k,v)

let cst v : linear = set Map.empty 0 v
let cst_int32 v : value = Semantics.Int(32,Range.point v),    cst v 
let cst_elt   v : value = Semantics.Elt(Some(Range.point v)), cst v 

let tryCst (a:linear) = 
  if   a.Count = 0                    then Some zero
  elif a.Count = 1 && a.ContainsKey 0 then Some a.[0] 
  else None

let linear_zero = cst zero
let linear_one  = cst one


let coerce_elt (s,a) = 
  match s with
  | Semantics.Int(_,r) -> Semantics.Elt(Some r),a 
  | _                  -> failwith "double coercion to element"

let neg ((s,a):value) = 
  Semantics.neg s, Map.map (fun i ai -> -ai) a

let add ((s0,a0):value) ((s1,a1):value) : value = 
  Semantics.(+) s0 s1,
  if a0.Count < a1.Count // trying not to be too quadratic
  then Map.fold (fun r i a0i -> set r i (get a1 i + a0i)) a1 a0 
  else Map.fold (fun r i a1i -> set r i (get a0 i + a1i)) a0 a1 

let sub ((s0,a0):value) ((s1,a1):value) : value =
  Semantics.(+) s0 (Semantics.neg s1),
  Map.fold (fun r i a1i -> set r i (get a0 i - a1i)) a0 a1

let addn ((s,a):value) n = 
  Semantics.shift s n, 
  set a 0 (get a 0 + n)

let subn v n = addn v -n

/// Multiply by a constant
let muln (a:linear) (n:integer) = 
  if n.IsZero  then linear_zero
  elif n.IsOne then a
  else Map.map (fun i ai -> ai * n) a       

let mulnV (s,a) n = Semantics.muln s n, muln a n

/// Add a wire c with equation a * b = c 
let mul' range (a: linear) (b: linear) : linear = 
  // should check range against semantics & insert truncation if required
  let r = newRoot() 
  let a,b = if a.Count < b.Count then b,a else a,b // V is more efficient than W
  Map.iter (setV r) a 
  Map.iter (setW r) b 
  let k = newWire { bank = Local; V = p0; W = p0; Y = p1 r }
  linear1 k

/// Multiply, trying to save the equation.
let mul ((s0,a0):value) ((s1,a1):value) : value = 
  let s = Semantics.(*) s0 s1
  s,
  match tryCst a0, tryCst a1 with 
  | Some m, Some n -> Log.record Log.Mlt0; cst (m*n)
  | None, Some n   -> Log.record Log.Mlt0; muln a0 n
  | Some m, None   -> Log.record Log.Mlt0; muln a1 m
  | _              -> Log.record Log.Mlt1; mul' (Semantics.range s)  a0 a1

/// Divide (assuming a field semantics; could be optimized) 
let div ((s0,a0):value) ((s1,a1):value) : value = 
  // add wire c with equation b * c = a 
  let r = newRoot() 
  Map.iter (setY r) a0 
  Map.iter (setV r) a1
  let k = newWire { bank = Local; V = p0; W = p1 r; Y = p0 }
  Semantics.elt, linear1 k


let linearBit = function
    | 0 -> linear_zero
    | 1 -> linear_one
    | _ -> failwith "not a bit"

let bit0 = Semantics.int1_zero, linear_zero
let bit1 = Semantics.int1_one , linear_one
let valueBit = function
    | 0 -> bit0
    | 1 -> bit1
    | _ -> failwith "not a bit"

// encoding of boolean operations
// we sometimes have to adjust the range to 0..1
let bitCoerce ((s,a) as v)=
  match s with 
  | Semantics.Int(1,r) -> Semantics.Int(1,Range.inter r Range.bit),a
  | _                  -> failwithf "not a bit: %A" v

let bitNeg v = sub (Semantics.int1, linear_one) v
let bitAnd a b = mul a b   
let bitOr  a b = bitCoerce (add a (sub b (bitAnd a b)))             // a + b - ab
let bitXor a b = bitCoerce (add a (sub b (mulnV (bitAnd a b) two))) // a + b - 2ab

(*
// other variants we considered
let bitAnd a b = 
  if   a = linear_zero || b = linear_zero then linear_zero
  elif a = linear_one then b
  elif b = linear_one then a 
  else mul a b
let bitOr  a b = bitNeg (bitAnd (bitNeg a) (bitNeg b))
*)

/// (a == 0) ? 1 : 0   
let zerop ((s,a) as v:value) = 
  // 4 different proof techniques
  let known, binary = 
    match s with 
    | Semantics.Elt None                         -> false, false
    | Semantics.Elt(Some r) | Semantics.Int(_,r) -> not(Range.within zero r), r = Range.bit
  Log.record(Log.Zerop (not(known || binary)))
  if known then bit0 
  elif binary then bitNeg v 
  else 
    // if a's semantics is IntN, we may need to truncate first,
    // in that case, we are left to prove that all low-order bits are 0
    // which we do recursively on the sum of bits
    // r0:  a * k0 = 1 - k1
    // r1:  a * k1 = 0          
    let r0 = newRoot()
    let r1 = newRoot()
    let k0 = newWire { bank = Local
                       V = p0; 
                       W = p1 r0; 
                       Y = p0  }
    let k1 = newWire { bank = Local
                       V = p0
                       W = p1 r1
                       Y = p0.Add(r0, minus1  ) }
    Map.iter (setV r0) a
    Map.iter (setV r1) a
    setY r0 0 one
    Semantics.int1, linear1 k1

/// prove an equality assertion with an equation    
let zeroAssert' (s,a:linear) =
    // equation: 0 * 0 = a
    let r = newRoot()
    Map.iter (setY r) a

(* a less efficient variant:
    // equation: a * 1 = 0
    let r = newRoot()
    Map.iter (setV r) a
    setW r 0 one
*)

/// prove an equality assertion by eliminating a wire
/// used e.g. to link an output wire to a computed result
/// (can rightfully fail at compile-time, e.g. if a is constant)
let zeroAssert (s,a:linear) =
    match Map.tryFindKey (fun k x -> wire.[k].bank = Local && (x = one || x = minus1)) a with
    | Some k when Lib.subst -> 
        // We first find a suitable intermediate wire in a, 
        // preferably with coefficient 1 or -1 and few non-zero roots
        let w = wire.[k]
        failwith "TODO:wire <- wire.Remove k"
        let b = a.Remove k
        let b = if a.[k] = one then Map.map (fun i xi -> -xi) b else b
        // we have c_k := b 
        // whenever c_k appeared in an equation, substitute it with b
        Map.iter 
          ( fun i xi -> 
              Map.iter (fun r vkr -> setV r i (vkr * xi)) w.V
              Map.iter (fun r wkr -> setW r i (wkr * xi)) w.W
              Map.iter (fun r ykr -> setY r i (ykr * xi)) w.Y )
          b
        // we will get "key k not found" if the variable is used elsewhere
        // we could keep a table recording c_k := b and substitute lazily 
        // alternatively, we could substitute only after generating the QAP.
    | _ -> zeroAssert' (s,a) 

/// turns a compile-time constant i into a wire k with equation 0 * 0 = a - i    
let endorse (i:integer) (rg:Range.r) = 
    let r = newRoot()
    let k = newWire { V = p0
                      W = p0
                      Y = p1 r 
                      bank = Local }
    setY r 0 -i
    Semantics.Int(32, rg), linear1 k

let bitAssert (a:linear) = 
    // equation: a * ( 1 - a ) = 0 
    let r = newRoot()
    Map.iter (fun i xi -> setV r i xi   ) a
    Map.iter (fun i xi -> setW r i (-xi)) a
    setW r 0 one

let newBit() = 
    let r = newRoot()
#if false 
    // equation: x * (1 - x) = 0  
    let k = newWire { bank = Local; range = range.bit  
                      V = p0.Add(r, minus1)
                      W = p1 r
                      Y = p0 }
    setV r 0 one
#else
    // equation: (x - 1) * x = 0  (tweaked to have less coefficients on W)
    let p1r = p1 r
    let k = newWire { V = p1r
                      W = p1r 
                      Y = p0 
                      bank = Local }
    setV r 0 minus1
#endif
    Semantics.int1, linear1 k 

(* old:
let bitwidth rg : int = // inefficient; will do for now.
    // we support two cases, could also do all-negatives 
    // (positive)       0 <= x < 2^n
    // (signed)  -2^(n-1) <= x < 2^(n-1)
    if rg.min < integer.Zero then
      max (log2 -rg.min) (log2 (max integer.Zero rg.max))  
    else
      log2 rg.max *)

let coerce_bit s ((s1,v):value) = Semantics.coerce_bit s s1, v
 
/// computes sum_(i<b) b[i]*2^i
let mergeBits' semantics (bs: value array) : value = 
    let b = bs.Length
    let mutable weight = one
    let mutable a = coerce_bit semantics bit0  
    for i = 0 to b-1 do
      a <- add a (mulnV (coerce_bit semantics bs.[i]) weight)
      weight <- weight * two
    a

/// generate split wires & equations based on range 
let splitVar' ((semantics,v) as a:value) : plan * value array =
    let (b,s,h) as p = planSplit (Semantics.range semantics)  
    p, match s with 
       | Unsigned -> 
          let a' : value = subn a (h <<< b)
          let bs : value[] = Array.init b (fun _ -> newBit()) // representing split b a'
          // equation: sum b[i]*2^i = a - h*2^b  
          let m = mergeBits' semantics bs
          zeroAssert (sub m a')

          let bits = 
            match semantics with 
            | Semantics.Int(n,_) -> // truncate or extend Int32s from b to n bits 
                if n <= b then Array.sub bs 0 n 
                else
                  let hs = Array.map valueBit (split (n - b) h)
                  Array.append bs hs

            | Semantics.Elt _    -> assert (h = zero); bs
                // let hs = Array.map valueBit (split (n - b) h) in Array.append bs hs

          if Lib.v  then 
            printfn "splitting %d bits, %A, unsigned %A" b semantics a 
            // printfn "splitting (%d bits, %A, unsigned) %A ---> %A" b semantics a bits
          bits

       | Signed -> 
          let a' = subn a ((two*h - one) <<< (b - 1))
          let bs = Array.init b (fun _ -> newBit()) // representing split b a'
          // equation: sum b[i]*2^i = a - h*2^b + 2^(b-1)
          zeroAssert (sub (mergeBits' semantics bs) a')
          let n = match semantics with | Semantics.Int(n,_) -> n 
                                       | Semantics.Elt _    -> printfn "signed split semantics?? %A : %A" a semantics; 32
          // truncate or sign-extend from b to n bits
          let c = bs.[b-1]
          let bits = 
            if n <= b then Array.sub bs 0 n
            else   
              let hs0 = split (n - b) (h - one) 
              let hs1 = split (n - b) h
              let hs = Array.map2 (fun h0 h1 -> add (valueBit h0) (mulnV c (bigint (h1 - h0)))) hs0 hs1
              Array.append bs hs 
          // flip the sign bit
          if b <= n then bits.[b-1] <- sub (valueBit 1) c 
          if Lib.v then 
            printfn "splitting (%d bits, signed): %A --%A--> %A" b a p bits
          bits   

let inline split_add (a: linear) (plan: plan,bs: value array) = 
    if splitCache.Count > 16 then splitCache.Clear();  // a crude way of preventing large caches
    splitCache.[a] <- (plan, bs) // save the cost of re-splitting
let inline split_lookup (a:linear) : plan * (value array) = splitCache.[a] :?> plan * (value array) 

let mergeBits s (bs: value array) : value =
    let a = mergeBits' s bs
    // split_cache a bs // no point recording, as we don't have the plan anyway
    a

let splitVar ((s,a:linear) as v) = 
    // have we already split this? 
    if splitCache.ContainsKey a then
      if Lib.v then printfn "re-using split %A" a
      let plan, bs = split_lookup (a:linear)
      Log.record(Log.Split(s, false, plan));
      bs
    else
      let plan, bs = splitVar' v
      split_add a (plan,bs);
      Log.record (Log.Split(s, true, plan));
      bs

let geq (s,_ as a) = 
    // the splitVar generates an hint
    // we could turn some into zerop
    let plan = planSplit (Semantics.range s)

    match plan, s with
    | (b,Signed,h), Semantics.Int(n,_) -> zerop (splitVar a).[n  - 1]  // the sign bit
    | (_,Unsigned,_), _                -> printf "warning: geq on unsigned value %A\n" s; 
                                          Log.record(Log.Split(s, false, plan)) // at proof-time, no need for a split.
                                          valueBit 1
    | _                                -> failwithf "geq on Elt? %A" a


let copyOutput ((s,a) as v:value) =
    // we split & truncate only when the value may exceed native int32 representations 
    match s with 
    | Semantics.Int(32,rg) -> 
        let required = rg.min <= -(one <<< 32) || (one <<< 32) <= rg.max  
        Log.record (Log.Truncate required)
        if required then 
          let r = mergeBits s (splitVar v)
          if Lib.v then printfn "truncating output %A to %A" a r
          r
        else v
    | _ -> failwith "only supporting int32 in this mode"

let wireOutput k ((s,a) as v: value) =
    // equation: a * 1 = x
    let r = newRoot()
    Map.iter (setV r) a
    setW r 0 one
    setY r k one 


///-------------------- experiments for automated QAP partitioning

(* TAKE 1: does not work because the resulting graph is too large 
   (with SHA1, each variable link 600-970 roots )

/// Print a weighed graph for gpmetis
let metis() =
  // we ignore the unit variable (will be replicated anyway)
  // we associate to each root a node, clustering all roots with equations touching each memory bank.
  
  // seeing variables as hyper-vertices, we care only about their sets of roots
  let roots (w:wire) : r Set = 
    let collect (a: r Set) r _ = a.Add r
    let a = Set.empty
    let a = Map.fold collect a w.V
    let a = Map.fold collect a w.W
    let a = Map.fold collect a w.Y
    //printfn "%A for %A" a w;
    a

  let nn = ref 0 
  let newNode() = incr nn; !nn
  let node : int[] = Array.zeroCreate nRoot // 0 = no node, then 1,2,... starting with banked roots.
  let vertices: Map<r*r,int> ref = ref Map.empty 
  let bankNode : (Map<string,int>) ref = ref Map.empty
  let vertex n0 n1 = // we insert the same vertex twice
    // printfn "(%3d, %3d)" n0 n1
    if not(n0 = n1) then 
      let v = match Map.tryFind (n0,n1) !vertices with | Some v -> v | None -> 0
      vertices := (!vertices).Add((n0,n1), v + 1)
  Map.iter
    (fun k (w:wire) -> 
      match w.bank with 
      | One         -> () // we ignore the unit variable (pervasive, will be replicated anyway)
      | Shared name -> 
          // we process these variables first; 
          // we bundle all roots using variables in the bank into a single node
          let n = 
            match (!bankNode).TryFind name with
            | Some n -> n
            | None -> let n = newNode() in bankNode := (!bankNode).Add(name,n); n 
          Set.iter
            (fun r -> let n' = node.[r] in if n' = n || n' = 0 then node.[r] <- n else failwithf "double-banked %d %s" r name)
            (roots w)
      | Local ->
          // the local variables are those generating vertices
          let rs = roots w
          let ns = Set.map (fun r -> (if node.[r] = 0 then node.[r] <- newNode()); node.[r]) rs
          printfn "%d with %d roots and %d nodes." k rs.Count ns.Count
          Set.iter (fun n0 -> Set.iter (vertex n0) ns ) ns )
    wire
  // we need to count the vertices
  let nv = (!vertices).Count / 2 
  printf "%d %d 011" !nn nv
  let _ = 
    Map.fold
      (fun a (n0,n1) v -> 
        if a <> n0 then printf "\n%d" (if n0 < 3 then 100 else 1)  
        printf " %d %d" n1 v
        n0)
      0
      !vertices
  printfn "\n -----------"      
*)

(* TAKE 2, using a hybrid graph *)

/// Print a weighed graph for gpmetis
let metis graphfile parts =
  // we ignore the unit variable (will be replicated anyway)
  // to each root, we associate a node, merging any two roots that share a bank variable.
  // to each local, we also associate a node
  
  // seeing variables as hyper-vertices, we care only about their sets of roots
  let roots (w:wire) : r Set = 
    let collect (a: r Set) r _ = a.Add r
    let a = Set.empty
    let a = Map.fold collect a w.V
    let a = Map.fold collect a w.W
    let a = Map.fold collect a w.Y
    //printfn "%A for %A" a w;
    a

  let nn = ref 0 
  let newNode() = incr nn; !nn
  let node : int[] = Array.zeroCreate nRoot // 0 = no node, then 1,2,... starting with banked roots.
  let bankNode : (Map<string,int>) ref = ref Map.empty
  let vertices: Set<int*int> ref = ref Set.empty 
  let weight = Array.zeroCreate (nRoot + nWire) // indexed by nodes, zero by default (variables)
  let vertex nk nr = // we insert the same vertex twice
    assert(weight.[nk] = 0 && weight.[nr] <> 0)
    vertices := (!vertices).Add(nk,nr).Add(nr,nk)
  wire_iter
    (fun k (w:wire) -> 
      match w.bank with 
      | One         -> () // we ignore the unit variable (pervasive, will be replicated anyway)
      | Shared name -> 
        
        // not sure this is a good idea after all...
        if true then
          // we process these variables first: 
          // we bundle all roots using variables in the bank into a single node
          
          let n = 
            match (!bankNode).TryFind name with
            | Some n -> n
            | None -> let n = newNode() in bankNode := (!bankNode).Add(name,n); n 
          Set.iter
            (fun r -> let n' = node.[r]
                      if n' = n || n' = 0 
                      then
                        if n' = 0 then weight.[n] <- weight.[n] + 1; // we are merging this root
                        node.[r] <- n 
                      else failwithf "double-banked %d %s" r name)
            (roots w)
      | Local ->
          // the local variables are those generating vertices
          let rs = roots w
          let ns = Set.map (fun r -> (if node.[r] = 0 then let n = newNode() in node.[r] <- n; weight.[n] <- weight.[n] + 1); node.[r]) rs
          let nk = newNode();
          // printfn "x%d --> %d has %d roots" k nk rs.Count
          Set.iter (vertex nk) ns 
          ) 
  // we need to count the vertices
  let nv = (!vertices).Count / 2

  let f = System.IO.File.CreateText graphfile
  fprintf f "%d %d 010" !nn nv
  let _ = 
    Set.fold
      (fun a (n0,n1) -> 
        if a <> n0 then fprintf f "\n%d" (weight.[n0] (* * 100 + 1 *) ) 
        fprintf f " %d" n1
        n0)
      0
      !vertices
  f.Close()

  // playing with graphviz
  let f = System.IO.File.CreateText (graphfile + ".gv")
  fprintf f "strict graph QAP {\n  graph [TopView=\"1\"];\n1 -- {" 
  let _ = 
    Set.fold
      (fun a (n0,n1) -> // we are emitting each label twice. 
        if a <> n0 
        then fprintf f " };\n%d -- { %d" n0 n1 
        else fprintf f ", %d" n1
        n0)
      0
      !vertices
  fprintf f "};\n}\n"
  f.Close()

  // we don't know how good the partition is until we read back the result
  // however, the hybrid-graph cut is a safe approximation.
  printfn "\n -----------"      
  // Array.iteri (fun r n -> printfn "node(%6d) = %6d, weight %5d" r n weight.[n]) node 

  // this can work only the second time around, after running: gpmetis -contig t.gp 2
  // beware: p is indexed from 0, not 1!
  
  let resultfile = graphfile + ".part." + string parts
  if System.IO.File.Exists resultfile 
  then
      let p = 
        try Array.map System.Int32.Parse (System.IO.File.ReadAllLines resultfile)
        with _ -> failwithf "Missing partition %s." resultfile
      
      let sizes = Array.zeroCreate parts 
      Array.iter(fun x -> sizes.[x] <- sizes.[x] + 1) p
      // check the solution...
      let n = Set.fold (fun a (n,m) -> if p.[n-1] = p.[m-1] then a else a + 1) 0 !vertices / 2
      printfn "Nodes partition %A, cutting %4d edges" sizes n
      let rsizes = Array.zeroCreate parts
      Array.iter (fun n -> rsizes.[p.[n-1]] <- rsizes.[p.[n-1]] + 1) node
      let shared w = let rs = roots w in let ps = Set.fold (fun s r -> Set.add p.[node.[r]-1] s) Set.empty rs in if ps.Count > 1 then Some rs else None 
      let nv = 
        wire_fold
          (fun a _k w -> 
            match w.bank, shared w with 
            | Local, Some rs -> a+1
            | _              -> a  )
          0
      printfn ""
      printfn "Roots %d --> %A, sharing %d variables" nRoot rsizes nv

  
///-------------------- printing and debugging support

let private factor x = if x = one then "" else string x
let private printTerm getxi = 
    printf "("
    if wire_fold 
            (fun first i w -> 
              match getxi w with 
              | None -> first 
              | Some (xi: bigint) -> 
                if i = 0       then printf "%s" (string xi) 
                elif first     then printf "%sx%d" (factor xi) i
                elif xi < zero then printf " - %sx%d" (factor (-xi)) i 
                else                printf " + %sx%d" (factor xi) i
                false) true  
    then printf "0"
    printf ")"

/// Dump the whole QAP
let print() =
  printf "--- The quadratic program has degree %d and size %d:\n" nRoot nWire  
  for r = 0 to nRoot - 1 do
    printf "  root %3d:  " r  
    printTerm (fun w -> w.V.TryFind r)
    printTerm (fun w -> w.W.TryFind r)
    printf " = "
    printTerm (fun w -> w.Y.TryFind r); 
    printf "\n" 
  done  

/// Check that all equations actually hold for those recorded values
let check (values: Field.t array) =
  for r = 0 to nRoot - 1 do
    let v = ref Field.t.zero
    let w = ref Field.t.zero
    let y = ref Field.t.zero
    wire_iter
      ( fun i t ->
          match t.V.TryFind r with 
          | Some a -> v := !v + Field.t.ofBig a * values.[i] 
          | None -> ()
          match t.W.TryFind r with 
          | Some a -> w := !w + Field.t.ofBig a * values.[i]
          | None -> ()
          match t.Y.TryFind r with 
          | Some a -> y := !y + Field.t.ofBig a * values.[i]
          | None -> ())
    if not (!v * !w = !y) then 
      printfn "broken equation at root %d" r
      printTerm (fun w -> w.V.TryFind r)
      printTerm (fun w -> w.W.TryFind r)
      printf " = "
      printTerm (fun w -> w.Y.TryFind r); 
      printf "\n" 
      Array.iteri
        ( fun i t ->
          match t.V.TryFind r with 
          | Some a -> printfn "x%i = %A\n" i values.[i] 
          | None -> ()
          match t.W.TryFind r with 
          | Some a -> printfn "x%i = %A\n" i values.[i] 
          | None -> ()
          match t.Y.TryFind r with 
          | Some a -> printfn "x%i = %A\n" i values.[i]
          | None -> ())
        wire


module Test =
   
  let split1 range = 
    let k = newWire { bank = Shared "input"; V = p0; W = p0; Y = p0 }
    let v = Semantics.Elt (Some range), linear1 k             
    let bits = splitVar v
    print()
    printfn "%A" bits
    bits 

  let split() = 
    let b1 = split1 { min = integer.Zero; max = bigint 7 }
    let b2 = split1 (Range.uint 32)
    ()