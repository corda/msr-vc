module Field

// This module imports EC arithmetic on unmanaged field elements.

open System.Runtime.InteropServices

// Using nativeint as a void pointer
// This might not be portable / could become deprecated later on!

[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private makeEncoding(int)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private getField(void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private getRandomFieldElt(void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private field_add(void *, void *, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private field_sub(void *, void *, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private field_mul(void *, void *, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private field_mul_alt (void *, void *, void *, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private field_div(void *, void *, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private field_exp(void *, int, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern int  private fieldEltEql(void * f1, void * f2, void * f)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern int private field_compare(void *, void *, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private makeFieldElt(int, void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private printFieldElt(void *, void*)
[<DllImport(@"FFLib.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)>] extern [<MarshalAs(UnmanagedType.LPStr)>] string private fieldEltToString(void *, void*)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private deleteFieldElt(void *, void*)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private getMod(uint64[], void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private getFull(void *, uint64[], void *)
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private setFull(void *, uint64[], void *)

/// various elliptic curve implementations
/// for convenience, we implicitly use a global, mutable "current encoding"
let ENCODING_DBG_PRIME = 1
let ENCODING_DBG_32    = 2
let ENCODING_ENIGMA_BN = 3
let ENCODING_ARITH_BN  = 4
let ENCODING_ARITH_CP  = 5
let ENCODING_ARITH_QAP  = 6
let ENCODING_ARITH_BN_DBG  = 7
let ENCODING_ARITH_CP_DBG  = 8
let encoding_help = 
    "    ENCODING_DBG_PRIME    = 1\n" +
    "    ENCODING_DBG_32       = 2\n" +
    "    ENCODING_ENIGMA_BN    = 3 (default)\n" +
    "    ENCODING_ARITH_BN     = 4\n" +
    "    ENCODING_ARITH_CP     = 5\n" + 
    "    ENCODING_ARITH_QAP    = 6\n" +
    "    ENCODING_ARITH_BN_DBG = 7\n"+
    "    ENCODING_ARITH_CP_DBG = 8\n"

/// Field elements fit in (at most) 4 uint64 words in F# 
let mutable nwords = 4
let mutable nbytes_per_elt:int64 = int64(8 * 4);
let mutable encoding: nativeint = System.IntPtr.Zero
let mutable field: nativeint    = System.IntPtr.Zero
let mutable p : bigint = bigint 0
let private word = Array.create nwords 0UL // buffer for all conversions

let two64 = bigint 1 <<< 64
let two32 = bigint 1 <<< 32 
let two31 = bigint 1 <<< 31 

let private wordBig() = 
  let mutable r = bigint word.[nwords - 1]
  for i = nwords - 2 downto 0  do
    r <- two64 * r + bigint word.[i]
  r

let hex (s:bigint) =
  // we normalize to nwords * 8 as we print mostly field elements; ugly.
  let target = nwords * 8
  let b = s.ToByteArray() // sometimes with a zero digit at the top
  let a : byte[] = Array.zeroCreate (max b.Length (nwords * 8))
  if b.[b.Length-1] = 0uy 
  then Array.blit b 0 a 0 (b.Length-1) 
  else Array.blit b 0 a 0 b.Length 
  Array.fold (fun s b -> sprintf "%02X%s" b s) "" a

let private ctr = ref 0 // for deterministic sampling of `random' elements

/// abstract field elements (the native implementation is private)
//[<CustomEquality; CustomComparison; StructuredFormatDisplay("{AsString}")>]
[<StructuredFormatDisplay("{AsString}")>]
type t(e:nativeint) = 

    // we will allocate a few constants once we have the encoding 
    static let notYet = t(System.IntPtr.Zero)
    static let mutable Zero:t  = notYet
    static let mutable One:t   = notYet
    static let mutable Two:t   = notYet
    static let mutable Two32:t = notYet

    let p = System.GC.AddMemoryPressure nbytes_per_elt; incr Lib.nf; e    
    member x.raw = p

    static member setConstants() = 
      Zero  <- t.ofInt32 0 
      One   <- t.ofInt32 1
      Two   <- t.ofInt32 2
      Two32 <- t.ofBig two32

    static member zero = Zero
    static member one = One
    static member two = Two
    static member two32 = Two32

    static member (+) (x:t, y:t) = t(field_add(x.raw,y.raw,field))
    static member (-) (x:t, y:t) = t(field_sub(x.raw,y.raw,field))
    static member (~-) x         = Zero - x  
    static member (*) (x:t, y:t) = t(field_mul(x.raw,y.raw,field))
    static member (/) (x:t, y:t) = if y = Zero then failwithf "divide by zero" else t(field_div(x.raw,y.raw,field))

    // this is *not* mod p, would benefit from FFlib support
    // used only for truncating elements to their low-order bits
    static member (%) (a:t,b:t)  = t.ofBig (a.Big % b.Big)

    static member expf (x:t) (i:int)= t(field_exp(x.raw, i, field))
    static member cmpf (x:t, y:t) = field_compare(x.raw,y.raw,field)
      
    static member eqlf (x:t, y:t)= fieldEltEql(x.raw,y.raw,field) <> 0

    static member random() = 
      if Lib.random 
      then t(getRandomFieldElt(field))
      else incr ctr; t.ofInt32 !ctr
    static member (<<<) (x,i)          = x * t.expf Two i
    static member ofInt64 (i : int64) : t = t.ofBig (bigint i) 
    static member ofInt32 (i : int32) : t = 
      // note that i is treated unsigned
      t(makeFieldElt(i,field))
    static member ofBig(n: bigint) : t =
      let mutable r = n % p
      if r.Sign = -1 then r <- p + r
      for i = 0 to nwords - 1 do 
        let j = p
        word.[i] <- uint64 (r % two64)
        r <- r >>> 64
      let x = t.ofInt32 0 // allocating
      setFull(x.raw,word,field)
      x

    member x.Big : bigint = 
      // getting an *unsigned* bigint out of a field element 
      getFull(x.raw,word,field)
      let i = wordBig();
      i 
      // consider getting signed one in some cases, e.g.
      // if i < p / bigint 2 then i else i - p
    
    member x.Bit i = 
      // inefficient... Get FFLib support? 
      getFull(x.raw,word,field); word.[i / 64] &&& (1UL <<< (i % 64)) > 0UL
    member x.Print() = printFieldElt(x.raw,field); printf " -- \n"
    override x.Equals o   = match o with :? t as y -> t.eqlf(x,y) | _ -> failwith "field equality"
    override x.GetHashCode() = failwith "no hash code"
    interface System.IComparable with
        member x.CompareTo y = match y with :? t as y -> t.cmpf(x,y) | _ -> failwith "field comparison"
    override x.Finalize() = System.GC.RemoveMemoryPressure(nbytes_per_elt); deleteFieldElt (x.raw,field)
    
    // a property that calls a method, used for %A printing of elements.
    member x.AsString = x.ToString()  
    override x.ToString() = fieldEltToString(x.raw, field)


/// sets the current field (implicitly used for all operations)
/// we need in particular to reset the field constants
let setEncoding e = 
  nwords   <- match e with
              | 1 | 2 | 7 | 8 -> 1 
              | 6             -> 2
              | 3 | 4 | 5     -> 4
              | _             -> failwith "invalid encoding" 
  nbytes_per_elt <- int64(8 * nwords)
  encoding <- makeEncoding e 
  field    <- getField encoding
  for i = 0 to nwords - 1 do word.[i] <- 0UL
  p <- (getMod(word,field); wordBig())
  t.setConstants()
  if Lib.v then
    printfn "setting encoding to %d; p = %s." e (hex p)
// initial value, moved to Main, to prevent a potential loading race
// setEncoding ENCODING_ENIGMA_BN 


module IO =

  // Reading and Writing binary files of integers and field elements, imported from C++ for interop

  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private file_open(void*, string, bool)
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern int private file_read_unsigned_char(void*, void*)
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern int private file_read_int(void*, void*)
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void private file_write_int(void*, void*, int)
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void* private file_read_field_elt(void*, void* )
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void private file_write_field_elt(void*, void*, void* )
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void private file_close(void*, void* )

  type file = H of nativeint

  // C++ flags
  let Read = true
  let Write = false

  let fopen name flag = 
    if Lib.v then printfn "opening %s %A" name flag
    H(file_open(encoding, Lib.path + name, flag))

  let read_int size (H fh): int =
    let i = 
      match size with
      |  8 -> file_read_unsigned_char(encoding, fh)  
      | 32 -> file_read_int(encoding, fh)  
      |  _ -> failwith "unsupported native read"
    // printfn "reading %2d-bit int %d" size i
    i

  let write_int (H fh) (i:int) = 
    // printfn "writing %d" i
    file_write_int(encoding, fh, i)  

  let elem (H fh): t = 
    let e = t(file_read_field_elt(encoding,fh))
    if not (e.Big < p) then failwith "file_read_field_elt: not in the field!" 
    e
  let write_elem (H fh) (e:t) = file_write_field_elt(encoding, fh, e.raw)

  let fclose (H fh): unit = file_close(encoding, fh)

  // for "plaintext" banks, assumed to be int32s.
  let save_array name a = 
    let h = fopen (name + ".data") Write
    if Lib.v then printfn "saving %A" a
    let down n : int = 
      if n > two31 then 
        let v = int (n - two32)
        if Lib.v then printfn "normalizing %A to %d" n v
        v
      else int n
    Array.iter (fun (e:t) -> write_int h (down e.Big)) a
    fclose h

  let load_array name size = 
    let h = fopen (name + ".data") Read
    let a = Array.init size (fun _ -> let i = read_int 32 h in t.ofInt32 i)
    fclose h
    a


/// basic sanity checks for FFLib interop
let test() = 
    printfn "p is %A" p
    printfn "1/16 = %A " (t.one / (t.ofInt32 16));
    printfn "1/40 = %A " (t.one / (t.ofInt32 40));

    let mutable x = t.ofInt64 12L
    let mutable y = bigint 12
    for i = 1 to 20 do
        printfn "%82s" (x.ToString())
        printfn "%82s" (y.ToString())
        x <- x * x
        y <- (y * y) % p
    done
    printfn "-1 is %A, 1/-1 is %A" (-t.one) (t.one / -t.one) 
    let y = t.one / x
    printfn "1/x is %A, their product is %A" y (x*y) 
    let z = ref x
    let sz = 16
    let n = 1 <<< sz
    Lib.time (sprintf "2^%d field multiplications" sz) 
        (fun () -> for i = 1 to n do z := !z * !z done) ()
    let r0 = !z
    z := x
    let z' = (!z).raw 
    Lib.time (sprintf "2^%d field multiplications (alt)" sz) 
        (fun () -> for i = 1 to n do field_mul_alt(z',z',z',field) done) ()
    let r1 = !z
    let x = ref x.Big
    Lib.time (sprintf "2^%d big multiplications mod p" sz)  
        (fun () -> for i = 1 to n do x := (!x * !x ) % p done) ()
    printfn "\n%s vs \n%s vs \n%s" (string r0) (string r1) (x.Value.ToString())
    let x  = t.ofInt64 34L
    let x' = t.ofInt64 33L + t.one
    printfn "(%s = %s) is %b" (string x) (string x') ((x:t) = x')