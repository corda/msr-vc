module Val

// This module defines two abstract domains for generating 
// and evaluating QAPs by symbolic interpretation 
// 
//     v is used for compiling
// Run.v is used for evaluating

type integer = bigint

open QAP 

module Native = // plain, fixed-sized, C integers

  [<StructuredFormatDisplay("{AsString}");StructuralEquality;StructuralComparison>]
  type a =
    | Local of int32           // local address on the data stack
    | Global of string * int32 // global address (infix) 
  with 
    static member (+) (a,j) =
      match a with 
      | Local i      -> Local(i+j)
      | Global (b,i) -> Global(b,i+j) 
    member a.AsString = 
      match a with 
      | Local i -> sprintf "$%d" i
      | Global (b,i) -> sprintf "%s+%d" b i

  type t = int32 // fixed, signed integers; we could also use int64
  let p = one <<< 32 
  let max = bigint t.MaxValue
  let zero = 0
  let one = 1
  
  /// converts from big integers to signed natives, with truncation.
  let truncate (n:integer) =
    let t = n % p 
    let r = int  (if t > max then t - p else t)
    if Lib.v && integer r <> n then printfn "Truncating %A to %d" n r
    r

  /// converts from signed native to positive big integer 
  /// (possibly towards bitwise operations)
  let integer (n:t) = 
    if n < 0 then bigint n + max else bigint n 

  /// bitwise n-bit decomposition.
  /// valid inputs range from -2^(n-1) to 2^n - 1
  let splitn zero one n (x:t)  =
    if n < 32 && (x < -(1 <<< (n-1)) || x > (1 <<< n)) then failwithf "%d-split overflow: %A" n x 
    let a = Array.zeroCreate n;
    // ugly, avoiding signed extension issues.
    let mutable x' = int64 x + (if x < 0 then 1L <<< 32 else 0L)
    for i = 0 to n-1 do
      a.[i] <- if x' % 2L = 0L then zero else one
      x'     <- x' >>> 1
    if x' <> 0L then printfn "leaving %A" x' 
    // it is Ok to leave bits out of negative numbers.
    if Lib.v then printfn "split-%d %A to %A" n x a
    a

(* some proposed code restructuring:
module NextTime =
  // the interpreter values, i.e. whatever goes into local registers
  type 'int value = 
    | Int  of 'int      // integer
    | Addr of Native.a  // pointer 
    | Void              // undefined and inaccessible

  // our abstract domain for integer values
  // we keep the type and range of variable integers, as truncations may be deferred
  // we don't need them for known integers, as only their operations are typed.
  // [not what we do so far e.g. (+) is overloaded]
  type semantics = 
    | I of int // i1 i8 i16 i32 i64
    | Q                     
  type number = 
    | Int of Native.t 
    | Elt of Field.t
  type repr = 
    | Arith of linear        // linear combination of QAP variables
    | Bits  of linear array  // vector of linear combinations (some of them 0 and 1)
    | Top                    // no QAP at top-level              
  type variable = semantics * QAP.range.t * repr
  type 'var symbolic =
    | Cte of number 
    | Var of 'var
  type vGen = (variable symbolic) value
  type vRun = (number symbolic) value
  // in particular, we can use smallish integers for variables at evaluation-time as long as there is no overflow 
  // how to represent banks? Just use their embedded representation, as we do for keys etc.
*) 


// We have two "abstract" representations of values, for 
// arithmetic and bitwise operations; we convert between them on demand.
// We cache those conversions, to share any bit decomposition
// (not clear how it will scale up)
//
// We also have concrete values: integers and pointers.
// We switch to large integers only for abstract arithmetic 
//
// We do some light dynamic typing, separating addresses from pointers
// Pointer arithmetic is limited to Pointer + Offset.
                
/// generation-time representation of values
[<StructuredFormatDisplay("{AsString}")>]
type v =  
  | Var  of value       // symbolic value (computed range)
  | Bits of QAP.Semantics.t * value array // vector of symbolic bits (range 0..1) 
  | Int  of Native.t    // concrete small integer
  | Elt  of Field.t     // concrete field element (aka large integer modulo q)
  | Address of Native.a // concrete data pointer, approximating runtime type safety
  | Void 
  with 
  override x.ToString() = 
    match x with 
    | Var (s,_)  -> if true then sprintf "? %A" s else "?"
    | Bits (s,bs)-> sprintf "%d bit %A" bs.Length s //"%A" bs
    | Int n      -> sprintf "%d" n
    | Elt e      -> sprintf "%A" e
    | Address p  -> sprintf "%A" p
    | Void       -> sprintf "void"
  member m.AsString = m.ToString()
  static member zero = Int Native.zero
  static member one  = Int Native.one
  static member two32 = Elt Field.t.two32

  member x.uint32 =  // ensures that the value is within 0.. 2^32 - 1
    match x.Arith with
    | Var a    -> x // we should disallow negative values here.
    | Elt e    -> Elt(e % Field.t.two32)
    | Int n    -> if n < 0 then Elt(Field.t.ofBig(Field.two32 + bigint n)) else x  // defensive... recheck CPP interface
    | _        -> failwithf "not an uint32: %A" x
  member x.native = 
    match x.Arith with 
    | Int n -> n
    | _        -> failwithf "not a concrete integer: %A" x

  /// must be a compile-time boolean
  member x.True = x.native <> Native.zero 

  member x.knownTrue  = x = v.one  // used by DataCache
  member x.knownFalse = x = v.zero // used by DataCache
  member x.unknown  = // unknown at compile-time?
    match x.Arith with Var _ | Bits _ -> true  | _ -> false

  member x.address = 
    match x with 
    | Address n -> n
    | _         -> failwithf "not a pointer: %A" x
  member x.big =
    match x with
    | Int n    -> integer n 
    | Elt e    -> e.Big
    |     _    -> failwithf "toVal %A" x  
  member x.elt = 
    match x with 
    | Int n    -> Field.t.ofInt32 n
    | Elt e    -> e
    | _        -> failwithf "not an element: %A" x
  member x.var =
    match x with 
    | Int n -> cst_int32 x.big
    | Elt e -> cst_elt   x.big
    | _     -> failwithf "x.var: x should be known, is %A" x
  member x.bits : value array =
    match x with 
    | Bits(s,bs) -> bs
    | Var a      -> splitVar a 
    | _          -> Native.splitn (valueBit 0) (valueBit 1) 32 x.native 
  member x.bit i : v =
    match x with 
    | Var ((s,_) as a)  -> Bits(s,(splitVar a)).bit i // trusting the cache; no need to split if we extend?
    | Bits(s,bs)        -> if i < bs.Length then Var (bs.[i]) else v.zero 
    | Int n when i > 31 -> v.zero // note that "1 >>> 32 = 1" ! MSDN says it is undefined if i > 31  
    | Int n             -> if (n >>> i) &&& 1 = 1 then v.one else v.zero
    | Elt e             -> if e.Bit i then v.one else v.zero
    | _                 -> failwithf "no bits in %A" x
   
  member x.Arith = // factoring out arithmetic conversions
    match x with 
    | Bits(s,bs) -> Var(mergeBits s bs)
    | x          -> x
  member x.isWire k = 
    match x with
    | Var (_,l) when l = linear1 k -> true
    | _  -> false
  static member (+) (x:v, y:v) : v = 
    match x.Arith, y.Arith with
    | Var a, Var b        -> Var(add a b)  
    | Var a, z | z, Var a -> Var(addn a z.big) // no overflow check yet
    | Address p, Int n 
    | Int n, Address p    -> Address (p + int n) 
    | Elt e, y            -> Elt(e + y.elt)
    | x, Elt e            -> Elt(x.elt + e)
    | Int n, Int m        -> Int(n + m)
    | _                   -> failwithf "%A + %A" x y 
  static member neg (x:v) = 
    match x.Arith with
    | Var a               -> Var(neg a)
    | Elt e               -> Elt -e
    | Int n               -> Int -n
    | _                   -> failwithf "neg not implemented: %A" x
  static member (-) (x:v, y:v) : v = 
    match x.Arith, y.Arith with 
    | Var a, Var b        -> Var(sub a b)  
    | Var a, z            -> Var(addn a (- z.big))
    | z, Var a            -> Var (sub z.var a) // Var(sub z.var a) // no overflow check yet
    | Address p, Int n    -> Address (p + (-n)) 
    | Elt e, y            -> Elt(e - y.elt)
    | x, Elt e            -> Elt(x.elt - e)
    | Int n, Int m        -> Int(n - m)

    // a hack, to tolerate testing null pointers
    | Void, Void          -> v.zero      
    | Address _, Void | Void, Address _ -> v.one
    | _                   -> failwithf "%A - %A" x y 
  static member (*) (x:v, y:v) : v = 
    match x.Arith, y.Arith with 
    | Var a, Var b        -> Var(mul a b)   
    | Var a, z | z, Var a -> Var(mulnV a z.big)
    | Elt e, y            -> Elt(e * y.elt)
    | x, Elt e            -> Elt(x.elt * e)
    | Int n, Int m        -> Int(n * m)
    | _                   -> failwithf "%A * %A" x y
  static member (/) (x:v, y:v) : v = // depends on the semantics! 
    match x.Arith, y.Arith with 
    | z, Var b            -> Var(div z.var b) 
    | Var a, z            -> Var(div a z.var)
    | Elt e, y            -> Elt(e / y.elt)
    | x, Elt e            -> Elt(x.elt / e)
    | Int n, Int m        -> Int(n / m)
    | _                   -> failwithf "%A / %A" x y

  static member (%) (x:v, y:v) : v = // only for compile-time integers for now 
    match x.Arith, y.Arith with 
    | Int n, Int m        -> Int(n % m)
    | _                   -> failwithf "%A / %A" x y
  
  static member bitwise (bit,bitwiseN,bitwiseB) (x,y) = 
    //if Lib.v then printfn "%A bitwise %A" x y 
    match x,y with 
    | Int n, Int m              -> Int (bitwiseN n m)
    | (Var(s,_) | Bits(s,_)), _
    | _, (Var(s,_) | Bits(s,_)) -> Bits(s, Lib.bits2 bit (valueBit 0) x.bits y.bits)
    | _                         -> failwithf "bitwise %A %A" x y

  static member And = v.bitwise (bitAnd,(&&&),(&&&)) 
  static member Or  = v.bitwise (bitOr,(|||),(|||))
  static member Xor = v.bitwise (bitXor,(^^^),(^^^))

  static member Shift (x:v, s:int) =  
   let r = 
    match x with 
    | Int n -> // we need an *unsigned* shift -- ugly. 
               Int(int32 (if s > 0 then uint32 n <<< s else uint32 n >>> -s))
    
    // disabled for SHA1
    // special case, undoing C optimization to save a bit split
    | Var a when 0 <= s && Lib.unshift -> Var(mulnV a (one <<< s))

    | Var (Semantics.Int(n,r) as sem,_) 
    | Bits(Semantics.Int(n,r) as sem,_)  -> 
               let bs = x.bits
               // shifting implicitly truncates to 32 bits.
               Bits(Semantics.uint32, // could be more precise
                 Array.init 32 
                   (fun i -> if 0 <= i - s && i - s < n then bs.[i-s] else valueBit 0)) 
   if Lib.v then printfn "%A <<< %d ---> %A" x s r
   r

  static member zerop (x:v) =
    match x.Arith with  
    | Var a    -> Var(zerop a)  
    | Elt e    -> if e = Field.t.zero then v.one else v.zero
    | Int n    -> if n =  Native.zero then v.one else v.zero 
    | _        -> failwithf "zerop %A" x
  static member gt0p (x:v) = v.ge0p (x - v.one)
  static member ge0p (x:v) =
    match x.Arith with 
    | Var a    -> Var(geq a)
    | Int n    -> if n >= 0 then v.one else v.zero
    | _        -> failwithf "no comparison on %A" x
  static member zeroAssert (x:v) = 
    match x.Arith with  
    | Var a    -> zeroAssert a
    | _        -> if x <> v.zero then failwithf "zeroAssert %A" x

  // specific to generation:
  static member endorse (i:integer) = Var (QAP.endorse i Range.uint32)
  static member coerce_elt (x:v) =
    match x with // there is no point coercing bits
    | Var a    -> Var(coerce_elt a) 
    | _        -> x
  static member nRoot() = 
    Log.record (Log.Nroot QAP.nRoot);
    if Lib.v then 
      printfn "%8d equations so far" QAP.nRoot
    Int(QAP.nRoot); // for debugging and adapting (for the latter we'll need a log entry)
  member x.bound r = 
    match x.Arith with 
    | Var (s,v) -> Var(Semantics.bound s r,v)
    | _           -> x
  member x.output =  
    match x.Arith with
    | Var a -> copyOutput a
    | _     -> x.var  
    // should we be more restrictive? 

let newWire b s = Var(s, linear1 (newWire { bank = b; V = p0; W = p0; Y = p0 }))

let ofWire b _ = newWire b Semantics.uint32

// these two are initialized at toplevel keygen
let mutable k0 = Map.empty
let mutable banks: (string * v[])[] = [||]
let bank name size = 
  k0 <- k0.Add(name,nWire)
  Array.init size (ofWire (Shared name))
let getbank name = snd (Array.find (fun(n,vs) -> n = name) banks)



module Run = 
/// evaluation-time representation of values; 
/// simpler, but we still need to separate Vars from the rest, 

 let mutable inter = Map.empty // should use an array instead
 let mutable nWire = 1 // not actually needed, but convenient for debugging.

 // emits a value towards crypto-proving; could be streamed.
 let log x = 
   // printfn "logging x%d = %A" nWire x
   inter <- inter.Add(nWire,x)
   nWire <- nWire + 1

 // The only change from the "Gen" value abstraction is that the inputs are known.
 // We represent them in the field, as they may temporary exceed native bounds 
 // We carefully log the intermediate values using side effects.
 type linear = Field.t
 let cst x   : linear = linear.ofBig x // was: ofInt64 (int64 x) // that makes a difference when x is negative. 
 let endorse i = let x = cst i in log x; x
 let add x y : linear = x + y
 let neg x   : linear = -x
 let sub x y : linear = x - y
 let mul (x:linear) (y:linear) : linear = 
   let z = x * y
   match Log.hint() with
   | Log.Mlt1 -> log z
   | Log.Mlt0 -> ()
   | h             -> failwithf "bad hint (expected Mlt) %A" h
   z

 let div (x:linear) (y:linear) : linear = 
   let z =
     if y = linear.zero 
     then printfn "div/0"; linear.zero // this can happen as we evaluate branches whose results are discarded
     else x / y
   log z; 
   z 

 let addn x n = x + (Field.t.ofBig n)
 let muln x n = x * (Field.t.ofBig n) 


 let bitAnd a b = mul a b 
 let bitOr  a b = add a (sub b (mul a b)) 
 let bitXor a b = add a (sub b (muln (mul a b) two))    
 let bitNeg a = sub linear.one a

 let zerop (x:linear) =
    match Log.hint() with
    | Log.Zerop false -> 
          if x   = linear.zero then linear.one
          elif x = linear.one then linear.zero
          else failwithf "not a bit: %A" x
    | Log.Zerop true ->
        if x = linear.zero then
          log(linear.zero)
          log(linear.one)
          linear.one
        else 
          log(linear.one / x)
          log(linear.zero)
          linear.zero
    | _ -> failwith "Bad Zerop hint"


 let zeroAssert (a:linear) = 
   if not (a = linear.zero) then failwithf "runtime assertion failed on %A" a

 let splitv n (x: integer)  =
  let a = Array.zeroCreate n;
  let mutable x = x
  for i = 0 to n-1 do
    a.[i] <- if x.IsEven then linear.zero else linear.one
    x     <- x >>> 1
  if Lib.v then printfn "Splitting-%d %A to %A" n x a
  a
  
 let linearBit = function
   | 0 -> linear.zero
   | 1 -> linear.one
   | _ -> failwith "not a bit"

 let private splitVar2 (a:linear) =
   match Log.hint() with 
   | Log.Split (semantics, fresh,((b,s,h) as plan)) -> 
       let bs, bits = QAP.splitInteger semantics plan a.Big 
       if fresh then 
         if Lib.v then 
           printfn "splitting (%d fresh bits): %A --%A--> %A" b a plan bits 
         Array.iter (fun x -> log (linearBit x)) bs
       plan, Array.map linearBit bits
    | h -> failwithf "bad hint (expected Split) %A" h

 let splitVar semantics (a:linear) = snd (splitVar2 a)

 let geq (a:linear) = 
   let plan, bits = splitVar2 a 
   // we rely on a representing a signed intN interger.
   match plan with
   | (_,Signed,_)   -> zerop bits.[31]
   | (_,Unsigned,_) -> (if a.Big > Field.two31 then failwithf "geq unsigned 0 <= %A" a); linear.one
   // printfn "0 <= %A aka %A is %b" a (a.Big - Field.p) (r = linear.one)

 /// unsigned merge (unclear when/how to get a smaller, signed merge)
 let mergeBits (bs: linear array) : linear =
    let n = bs.Length
    let mutable a = bs.[n-1] 
    for i = n-2 downto 0 do 
      a <- add (muln a two) bs.[i] 
    a

 let copyOutput (a: linear) : linear = 
    match Log.hint() with
    | Log.Truncate true -> let a' = mergeBits (splitVar Semantics.int32 a)
                           if Lib.v then printfn "truncated %A %A" a a'
                           a'
    | Log.Truncate _    -> a
    | h                 -> failwith "bad hint (expected Truncate) %A" h


// the rest is copied verbatim; could also use an abstract interface (slower?)
 [<StructuredFormatDisplay("{AsString}")>]
 type v =  
  | Var  of linear       // symbolic value (computed range)
  | Bits of linear array // vector or symbolic bits (range 0..1) 
  | Int  of Native.t     // concrete small integer
  | Elt  of Field.t      // concrete field element (aka large integer modulo q)
  | Address of Native.a  // concrete data pointer, approximating runtime type safety
  | Void 
  with 
  override x.ToString() = 
    match x with 
    | Var a     -> if a.Big < Field.two32 then sprintf "%A" a.Big 
                   // if a.Big - Field.p > -Field.two32 then sprintf "%A" (a.Big - Field.p)
                   else sprintf "%A" a
    | Bits bs   -> sprintf "%d bits" bs.Length //"%A" bs
    | Int n     -> sprintf "%d" n
    | Elt e     -> sprintf "%A" e
    | Address p -> sprintf "$%A" p
    | Void      -> sprintf "void"
  member m.AsString = m.ToString()
  static member zero = Int Native.zero
  static member one  = Int Native.one
  static member two32 = Elt Field.t.two32

  member x.isBit = 
    let r = match x.Arith with
            | Var e 
            | Elt e -> let b = e.Big in bigint 0 <= b && b <= QAP.one
            | Int n -> 0 <= n && n <= 1
            | _ -> false
    match Log.hint() with
    | Log.IsBit false       -> false
    | Log.IsBit true when r -> true
    | h -> failwithf "Was expecting a bit, but isn't %A" h

  member x.uint32 =  // ensures that the value is within 0.. 2^32 - 1
    match x.Arith with
    | Var a    -> assert (a.Big < Field.two32); x 
    | Elt e    -> Elt(e % Field.t.two32)
    | Int n    -> if n < 0 then Elt(Field.t.ofBig(Field.two32 + bigint n)) else x  // defensive... recheck CPP interface
    | _        -> failwithf "not an uint32: %A" x

  member x.native = 
    match x with 
    | Int n -> n
    | Elt e -> printfn "casting int %A" x; int e.Big
    | _        -> failwithf "not a concrete integer: %A" x
  member x.True = x.native <> Native.zero 
  member x.knownTrue  = match x.Arith with Int n -> n <> 0 | _ -> false
  member x.knownFalse = match x.Arith with Int n -> n = 0  | _ -> false
  member x.unknown    = match x.Arith with Int _ | Elt _ -> false  | _ -> true
  member x.address = 
    match x with 
    | Address n -> n
    | _         -> failwithf "not a pointer: %A" x
  member x.big =
    match x with
    | Int n    -> integer n 
    | Elt e    -> e.Big
    |     _    -> failwithf "toVal %A" x  
  member x.elt = 
    match x with 
    | Int n    -> Field.t.ofInt32 n
    | Elt e    -> e
    | _        -> failwithf "not an element: %A" x
  member x.var = cst x.big
  member x.bits =
    match x with 
    | Bits bs  -> bs
    | Var a    -> splitVar Semantics.elt a 
    | _        -> Native.splitn linear.zero linear.one 32 x.native 
  member x.bit i : v =
    match x with 
    | Var a    -> (Bits (splitVar Semantics.elt a)).bit i // trusting the cache; no need to split if we extend?
    | Bits bs  -> if i < bs.Length then Var (bs.[i]) else v.zero 
    | Int n when i > 31 -> v.zero
    | Int n    -> Int((n >>> i) &&& 1)  // dodgy for negative numbers
    | Elt e    -> Int(if e.Bit i then 1 else 0)
    | _        -> failwithf "no bits in %A" x

  member x.Arith = // factoring out arithmetic conversions
    match x with 
    | Bits bs -> Var(mergeBits bs)
    | x       -> x
  static member (+) (x:v, y:v) : v = 
    match x.Arith, y.Arith with
    | Var a, Var b        -> Var(add a b)  
    | Var a, z | z, Var a -> Var(addn a z.big) // no overflow check yet
    | Address p, Int n 
    | Int n, Address p    -> Address (p + int n) 
    | Elt e, y            -> Elt(e + y.elt)
    | x, Elt e            -> Elt(x.elt + e)
    | Int n, Int m        -> Int(n + m)
    | _                   -> failwithf "%A + %A" x y 
  static member neg (x:v) = 
    match x.Arith with
    | Var a               -> Var(neg a)
    | Elt e               -> Elt -e
    | Int n               -> Int -n
    | _                   -> failwithf "neg not implemented: %A" x
  static member (-) (x:v, y:v) : v = 
    match x.Arith, y.Arith with 
    | Var a, Var b        -> Var(sub a b)  
    | Var a, z            -> Var(addn a (- z.big))
    | z, Var a            -> Var(sub z.var a) // no overflow check yet
    | Address p, Int n    -> Address (p + -n) 
    | Elt e, y            -> Elt(e - y.elt)
    | x, Elt e            -> Elt(x.elt - e)
    | Int n, Int m        -> Int(n - m)

        // a hack, to handle testing null pointers 
    | Void, Void          -> v.zero
    | Address _, Void | Void, Address _ -> v.one

    | _                   -> failwithf "%A + %A" x y 
  static member (*) (x:v, y:v) : v = 
    match x.Arith, y.Arith with 
    | Var a, Var b        -> Var(mul a b)   
    | Var a, z | z, Var a -> Var(muln a z.big)
    | Elt e, y            -> Elt(e * y.elt)
    | x, Elt e            -> Elt(x.elt * e)
    | Int n, Int m        -> Int(n * m)
    | _                   -> failwithf "%A + %A" x y 
  static member (/) (x:v, y:v) : v = 
    match x.Arith, y.Arith with 
    | z, Var b            -> Var(div z.var b) 
    | Var a, z            -> Var(div a z.var)
    | Elt e, y            -> Elt(e / y.elt)
    | x, Elt e            -> Elt(x.elt / e)
    | Int n, Int m        -> Int(n / m)
    | _                   -> failwithf "%A / %A" x y

  static member (%) (x:v, y:v) : v = // only for compile-time integers for now 
    match x.Arith, y.Arith with 
    | Int n, Int m        -> Int(n % m)
    | _                   -> failwithf "%A / %A" x y

  static member bitwise (bit,bitwiseN,bitwiseB) (x,y) = 
    match x,y with 
    | Int n, Int m        -> Int (bitwiseN n m)
    | _                   -> Bits (Lib.bits2 bit linear.zero x.bits y.bits)
  static member And = v.bitwise (bitAnd,(&&&),(&&&)) 
  static member Or  = v.bitwise (bitOr,(|||),(|||))
  static member Xor = v.bitwise (bitXor,(^^^),(^^^))
    
  static member Shift (x:v, s:int) = 
   let r = 
    match x with 
    | Int n -> Int(int32(if s > 0 then uint32 n <<< s else uint32 n >>> -s))
    
    // disabled for SHA(as this grows the representation by 30 bits), reconsider.
    | Var a when 0 <= s && Lib.unshift -> Var(muln a (one <<< s))

    | _     -> let bs = x.bits
               Bits(
                 Array.init 32 
                   (fun i -> if 0 <= i - s && i - s < min bs.Length 32 then bs.[i-s] else linear.zero)) 
   if Lib.v then printfn "%A <<< %d ---> %A" x s r
   r

  static member zerop (x:v) =
    match x.Arith with  
    | Var a    -> Var(zerop a)  
    | Elt e    -> if e = Field.t.zero then v.one else v.zero
    | Int n    -> if n =  Native.zero then v.one else v.zero 
    | _        -> failwithf "zerop %A" x
  static member gt0p (x:v) = v.ge0p (x - v.one)
  static member ge0p (x:v) =
    match x.Arith with 
    | Var a    -> Var(geq a)
    | Int n    -> if n >= 0 then v.one else v.zero
    | _        -> failwithf "no comparison on %A" x
  static member zeroAssert (x:v) = 
    match x.Arith with  
    | Var a    -> zeroAssert a
    | _        -> if x <> v.zero then failwithf "zeroAssert %A" x

  // specific to evaluation:

  static member endorse x = Var(endorse x)
  static member coerce s (x:v) = x
  static member nRoot() = 
    match Log.hint() with 
    | Log.Nroot n -> Int n
    | _           -> failwith "Bad nRoot hint"
  member x.bound rg =
    match x.Arith with 
    | Var a -> if Lib.v && not(QAP.Range.within a.Big rg) then printfn "bound Var: %A not in %A" x rg
    | Int n -> if Lib.v && not(QAP.Range.within (bigint n) rg) then printfn "bound Int: %A not in %A" x rg
    | Elt e -> printfn "bound Elt %A in %A" e rg
    | _     -> printfn "bound ??? %A in %A" x rg
    x
  static member witness x = 
    match x with 
    | Var w -> w
    | _     -> failwithf "no witness in %A" x
  member x.output: linear = 
    match x.Arith with
    | Var n -> copyOutput n
    | _     -> x.var


 let enter() =
   // not great
   inter <- Map.empty
   nWire <- 1 
  
 let newIOs isize osize = 
   nWire <- 1 + isize