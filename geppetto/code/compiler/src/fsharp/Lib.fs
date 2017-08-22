module Lib

/// Kitchen sink, with global flags and utility functions.
 

let mutable t  = true      // print timings
let mutable bv = true      // verbosity of scheduling (~ VERBOSE in geppetto.h)
let mutable v  = false     // verbosity at top-level
let mutable vv = false     // higher verbosity (prints the mQAP, one coeff at a time)
let mutable i  = false     // verbosity in the interpreters 
let mutable e  = false     // print local QAPs in keygen
let mutable h  = false     // print histogram of variable widths at compile- and run-time
let mutable check = false  // check all local QAP equations hold before proving (slow)
let mutable random = false // does elem_rand sample randomly or from a fixed seed? 
let mutable fork = true    // do we compile branching on runtime values?
let mutable subst = false  // eliminate equated variables (experimental)
let mutable raw = false    // print timings in a machine-readable format
let mutable w = false      // whole-program mode (otherwise use command-line scripting) 
let mutable gc = false     // aggressive preemptive GCs to stabilize experimental results

let mutable iloop = 0 // loop counter for the "inner" verifiable code (deprecated)
let nf = ref 0        // number of F#-allocated field elements

// turn << s into * 2^s on variables; this undoes an LLVM optimization to avoid a bit decomposition
// however, it makes later bit decomposition more expensive by s bits. 
let mutable unshift = false 

let mutable expect_to_fail = true // flag passed to FFLib for crypto checks

//---- files 

/// there are two different binary formats for keys, proofs, and commitments. 
/// from F#, we can choose the format when saving, but we can only read back the rich format
/// from interpreted C, we can only read the simple format, e..g for bootstrapping
type fileFormat = Simple | Rich
let mutable path = "./build/"
let mutable save : fileFormat option = None  // save proof elements; Simple towards bootstrapping
let mutable emit : string option = None      // print C definitions at keygen time  


module CPP = // unmanaged configuration & final cleanup

  open System.Runtime.InteropServices
  [<DllImport(@"FFLib.dll")>] extern void private configure(int, int, int, int)
  [<DllImport(@"FFLib.dll")>] extern void cleanup()
  [<DllImport(@"FFLib.dll")>] extern int64 cGetRDTSC()
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void  private print_current_field_stats ()
  [<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern [<MarshalAs(UnmanagedType.BStr)>] string field_stats();

  // defaults copied from Config.cpp
  let mutable v     = 1 // verbosity for C++
  let mutable maxGB = 6 // size for exponentiation tables, in GBs.
  let mutable fpc   = 0 // pre-allocates bigger exponentiation tables 
 
  let set() = configure(maxGB,fpc, v, if raw then 1 else 0)

  

open System.Threading
       
let private time0 = System.DateTime.Now

/// functional wrapper for measuring time and other statistics 
let time (s:string) f x = 
  if t then
    //if raw then
    //  printf "%s: " (s.Replace(' ', '_'))
    //else 
      //printf "%30s... " s
    if not raw && v then
      printfn "%30s... " s    
    if gc then
      //Each of the following lines fails to produce stable timing results:
      // 1) System.GC.Collect()
      // 2) for i = 0 to System.GC.MaxGeneration do
      //      System.GC.Collect(i, System.GCCollectionMode.Forced, true)
      // 3) let mem = System.GC.GetTotalMemory true
      // 4) ignore (System.Runtime.GCSettings.LatencyMode = System.Runtime.GCLatencyMode.SustainedLowLatency)
      // This pair works
      System.GC.Collect()
      Thread.Sleep(1000)    

    //let fs = CPP.field_stats()
    ignore stdout.Flush  // try to avoid crossing the .NET and C++ output streams

    let collections0 = [for i in 0..System.GC.MaxGeneration ->  System.GC.CollectionCount i]
    let proc = System.Diagnostics.Process.GetCurrentProcess()
    nf := 0;
    let processor_start = proc.TotalProcessorTime
    let timer = new System.Diagnostics.Stopwatch() // now use a fresh one, to support nested timings
    timer.Start()
    let cycle_start = CPP.cGetRDTSC()
  
    let r = f x
  
    let cycles         = CPP.cGetRDTSC() - cycle_start
    timer.Stop()
    let processor_time = proc.TotalProcessorTime - processor_start
    let collections    = List.mapi (fun i start -> System.GC.CollectionCount i - start) collections0
    if raw then
      printfn "%s: %d %d %d %d %A" (s.Replace(' ', '_')) timer.ElapsedMilliseconds cycles timer.ElapsedTicks processor_time.Ticks collections
    else
      // simplified, as most numbers were redundant      
      // this is forcing a full GC, may affect performance      
      let m = string (System.GC.GetTotalMemory gc / 1000000L)
      // printfn "%25s %10d mS%10d mS%18d cycles   %13d ticks   %s%11d elts%6s MB after %s GCs" 
      printfn "%35s %10d mS%18d cycles     %s%6s   MB after %s GCs" 
        s
        timer.ElapsedMilliseconds 
        //processor_time.Milliseconds
        cycles
        //timer.ElapsedTicks 
        ((System.DateTime.Now - time0).ToString "hh':'mm':'ss" )
        //!nf
        m
        (String.concat "/" (List.map string collections))    
    //let fs = CPP.field_stats()
    r
  else f x

let pause x = ignore(System.Console.ReadLine()); x

let rec find k = function 
  | (x,v)::_ when x = k -> v
  | _::m -> find k m
  | []   -> failwithf "%A not found" k

let extend v0 (a: 'a array) i =
  if i < a.Length then a.[i] else v0

let bits2 f v0 (a: 'a array) (b: 'a array) = 
  // this late truncation assumes bits2 is used only on int32s
  let s = min 32 (max a.Length b.Length) 
  Array.init s
    (fun i -> 
      f (if i < a.Length then a.[i] else v0)
        (if i < b.Length then b.[i] else v0))

// In case of common keys: we take values from m2
let mapUnion m1 m2 = Map.fold (fun m k v -> Map.add k v m) m1 m2

let resize a n = 
  let r = ref a 
  array.Resize(r,n)
  !r 

open System.IO
open System.Runtime.Serialization.Formatters.Binary
open System.Runtime.Serialization
open System.Xml.Serialization

let xml = false
let toFile (name:string) x = 
  let name = path + name + ".mqap"
  if xml then    
    let fs = new StreamWriter(name)
    let t = x.GetType()
    let format = new XmlSerializer(t)
    format.Serialize(fs,x)
    fs.Close()
  else
    let fs = new FileStream(name, FileMode.Create)
    let format = new BinaryFormatter()
    format.Serialize(fs,x)
    fs.Close()

let fromFile name = 
  let name = path + name + ".mqap"
  let fs = new FileStream(name, FileMode.Open)
  let format = new BinaryFormatter()
  let v = format.Deserialize(fs)
  fs.Close()
  v
    
  
(*
/// Colored printf
let cprintf c fmt = 
    Printf.kprintf 
        (fun s -> 
            let old = System.Console.ForegroundColor 
            try 
              System.Console.ForegroundColor <- c;
              System.Console.Write s
            finally
              System.Console.ForegroundColor <- old) 
        fmt
*)
