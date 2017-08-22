module Main

// Geppetto supports two separate programming styles: 
// * CLI-based, with a sequence of detailed actions 
// * program-based, with schedule & details collected from a whole C program

open System.IO

/// ad hoc parsing for command-line arguments: num[,num]* with no space
let randomInputs sz = 
    let rnd = System.Random(42)
    let sample _ =  rnd.Next(0, System.Int32.MaxValue) // (System.Int32.MinValue, System.Int32.MaxValue)
    Array.init sz sample (* % (1 <<< 20) *) 

let parseArgs n (args: string list) =
  let subparse (s:string) = 
    if s.StartsWith "random*" 
    then randomInputs (System.Int32.Parse (s.Substring 7))
    else Array.ofList [System.Int32.Parse s]
  match args with 
  | seq::_ when n = 0 && seq.StartsWith "-" -> [||], args
  | []     when n = 0                       -> [||], [] 
  | "random"::args -> randomInputs n, args
  | seq::args -> 
        let vs = 
          try Array.collect subparse (seq.Split ',')
          with _ -> failwithf "bad arguments: %s" seq 
        if vs.Length <> n then failwithf "expected %d arguments, got [| %s |]" n seq
        vs, args
  | _ -> failwithf "not enough integer arguments: %A" args

let parseStrings (args: string list) =
  match args with 
  | seq::_ when seq.StartsWith "-" -> [||], args
  | []                             -> [||], [] 
  | seq::args -> seq.Split ',', args

open System.Runtime.InteropServices
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern void enc_count_reset()
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern int enc_count_total_get()
[<DllImport(@"FFLib.dll", CallingConvention = CallingConvention.Cdecl)>] extern int enc_count_nonzero_get()

open Bank // defining our global state 
   
/// compile an outsourced function
let compile args name = 
  let resultType, formals, body = VM.Code.findFunDef ("@"+name)

  // build compile-time actual arguments, either compile-time constants or banks
  let actuals, // actual compile-time arguments 
      cfgs,    // including those compile-time values from the command line
      args,    
      banks    // ...and these banks
    = 
    List.fold 
      (fun (actuals,cfgs,args,banks) (t,name) -> 
        match t with 
        | LLVM.Ptr(LLVM.Tvar n) when n.StartsWith Bank.prefix 
              -> let b = bank (n.Substring (Bank.prefix.Length))
                 actuals @ [name,Bank.Bank b], cfgs, args, banks @ [b]
        | LLVM.Ptr t 
              -> // we need concrete values for t, for now only on the command line 
                 let size = VM.Code.tSize t 
                 let vals, args = parseArgs size args
                 actuals @ [name,Bank.Cfg(vals)], cfgs @ [vals], args, banks
        | _ -> failwithf "argument type %A" t
      )
      ([],[],args,[])
      formals  

  // we would need something more general for functions returning multiple banks.
  let banks, result = 
    match resultType with  
        | LLVM.Ptr(LLVM.Tvar n) when n.StartsWith Bank.prefix 
              -> let b = bank (n.Substring (Bank.prefix.Length))
                 banks @ [b], Some b
        | _   -> banks, None

  let banks = Array.ofList banks
  System.Array.Sort(banks) 
  let kb,_ = Array.fold 
               (fun (a:Map<string,int>,k) (b:b) -> a.Add(b.name,k), k + b.spec.Length) 
               (Map.empty,1) banks
  Val.k0 <- kb

  if Lib.v then printfn "Compiling %s with compile-time arguments %A and banks %A" name cfgs kb  
  let qap = 
    Lib.time 
      "QAP generation" 
      (QAPGen.eval ("@"+name) actuals) banks
  
  let f = { name = name; banks = banks; kb = kb; args = cfgs; qap = qap; k0 = -1 ; c0 = Proof.zeros; proof = Map.empty }
  s.functions <- s.functions @ [f]
  // if Lib.v then 
  let (nv,nw,ny) = QAP.nCoeff qap
  printfn "compiled %s (%d steps) to %d equations (%d V + %d W + %d Y coefficients) on %d variables from banks %A and %d locals" name !QAPGen.nop qap.nRoot nv nw ny qap.nWire f.kb (rho f)
  QAPGen.nop := 0; 
  
  // stats on the range of internal values (informing exponentiation tables)
  if Lib.h then 
    let histogram  = Array.zeroCreate 255; 
    QAP.wire_iter 
      (fun _ (w: QAP.wire) -> 
        if w.bank = QAP.Local then 
          let r = failwith "disabled" // was: w.range.max - w.range.min
          let i = min 253 (if r = QAP.zero then 0 else (1 + QAP.log2 r)) 
          histogram.[i] <- histogram.[i] + 1) 
    printfn "#bit #variable\n"
    Array.iteri (fun i n -> if n <> 0 then printfn "%3d %7d" i n) histogram

  args


let rec command (args: string list) = 
  match args with 
  | [] -> (if Lib.v then printfn "Done."); 0
  | "-r"::args      -> Lib.raw <- true; Lib.CPP.set(); command args
  | "-t"::args      -> Lib.i  <- true; command args 
  | "-v"::args      -> Lib.v  <- true; command args 
  | "-vv"::args     -> Lib.vv <- true; command args 
  | "-e"::args      -> Lib.e  <- true; command args 
  | "-h"::args      -> Lib.h  <- true; command args 
  | "-w"::args      -> Lib.w  <- true; command args 
  | "--gc"::args    -> Lib.gc <-true; command args
  | "--fflib"::"GB" ::n::args -> Lib.CPP.maxGB <- int n; Lib.CPP.set(); command args 
  | "--fflib"::"v"  ::n::args -> Lib.CPP.v     <- int n; Lib.CPP.set(); command args 
  | "--fflib"::"fpc"::n::args -> Lib.CPP.fpc   <- int n; Lib.CPP.set(); command args 
  | "-q"::args      -> Lib.t <- false; Lib.v <- false; Lib.i <- false; Lib.e <- false; Lib.bv <- false;
                       Lib.CPP.v <- 0; Lib.CPP.set(); command args 
  | "--no-fork"::args -> Lib.fork <- false; command args
  | "--format"::"Simple"::args -> Lib.save <- Some Lib.Simple; command args
  | "--format"::"Rich"::args   -> Lib.save <- Some Lib.Rich; command args
  | "--check"::args -> Lib.check <- true; command args
  | "--encoding"::n::args -> Field.setEncoding (int n); command args
  | "--opt"::"unshift"::args -> Lib.unshift <- true; command args

  | "--test"::"field"   ::args -> Field.test(); command args
  | "--test"::"split"   ::args -> QAP.Test.split(); command args 

  | "--test"::"keygen"  ::args -> Proof.Test.basic(); command args
  | "--test"::"simple"  ::args -> Proof.Test.simple(); command args
  | "--test"::"fg"      ::args -> Proof.Test.fg(); command args
  //  | "--test"::"add"     ::args -> Proof.Test.add(); command args
  //  | "--test"::"layered" ::args -> Proof.Test.layered(); command args
  //  | "--test"::"poly"    ::args -> Proof.Test.poly(); command args

  | "--test"::"cfgllvm" ::args -> TestCfgAnnotator.llvmPrintIfElse (); command args
  | "--test"::"cfg1"    ::args -> TestCfgAnnotator.handWrittenTests(); command args
  | "--test"::"cfg2"    ::args -> match parseArgs 2 args with [|gs;gc|],args -> TestCfgAnnotator.randomTests gs gc; command args
                                                            | _ -> failwith "Parse error"

  | "--path"::path::args -> Lib.path <- path; command args                                   
  | "--emit"::name::args -> Lib.emit <- Some name; command args

  | "--load"::file::args -> // loading LLVM assembly code
        if VM.code = [] then 
          VM.Code.load file
          if Lib.v then 
            printfn "Loaded %d definitions from %s" VM.code.Length file
          command args
        else failwithf "Multiple LLVM inputs are not supported: %s" file

  | "--load-qap"::file::datafile::args -> // loading, proving, verifying some external QAP (for PinoccioQ) 
        Pinoccioq.run file datafile
        command args

  | "--save"::mode::filename::args -> // saving code & keys, between compilation and execution

        // annoyingly, the .NET binary serializer fails on some Maps
        // so we need some pre-processing; we may also serialize in smaller chunks

        // there is no need to save the QAP for the verifier.
        let bformat0 (s:string) (b:b) = Map.toArray b.task
        let bformat1 (s:string) (b:b) = { b with task = Map.empty }
        let format f : string * int array list * string array                                  * QAP.f option                                            * int = 
                       f.name ,         f.args , Array.map (fun (b: Bank.b) -> b.name) f.banks , (if mode="verifier" then None else Some(QAP.format f.qap)), f.k0
        let old = VM.code, Map.map bformat0 s.banks, Map.map bformat1 s.banks, List.map format s.functions, s.banksize  
        // printfn "%A" old
        Lib.toFile filename old
        if mode = "prover"   || mode = "both" then Proof.Prove.saveKey  Lib.Rich filename s.ek.Value
        if mode = "verifier" || mode = "both" then Proof.Verify.saveKey Lib.Rich filename s.vk.Value
        command args  

  | "--resume"::mode::filename::args -> // resuming code & keys
        if mode = "prover"   || mode = "both" then s.ek <- Some(Proof.Prove.loadKey filename)
        if mode = "verifier" || mode = "both" then s.vk <- Some(Proof.Verify.loadKey filename)
        let bformat1 (b0 : Map<string,(int * verifierTask)[]>) (s:string) (b:b) = { b with task = Map.ofArray b0.[s] }
        let format (n,a,b,o,m) = 
          { name =n
            args = a
            banks = Array.map (fun x -> s.banks.[x]) b
            qap = match o with Some q -> QAP.parse q | None -> QAP.dummy
            k0 = m
            kb = Map.empty
            c0 = Proof.zeros
            proof = Map.empty }
        let code, banks0, banks1, functions, banksize = Lib.fromFile filename :?> LLVM.decl list * Map<string,(int * verifierTask)[]> * Map<string,b> * (string * int array list * string array * QAP.f option * int) list * int
        if VM.code = [] then VM.code <- code else failwith "Multiple LLVM inputs are not supported"
        s.banks     <- Map.map (bformat1 banks0) banks1
        s.banksize  <- banksize
        s.functions <- List.map format functions
        command args  

  | "--input"::name::"from"::filename::args -> 
        let b = bank name 
        let t = LLVM.Tvar("%struct." + "bank_" + name)
        let h = Field.IO.fopen filename Field.IO.Read
        let pos = ref 0
        let readI sz = 
          let i = Field.IO.read_int sz h
          b.value.[!pos] <- Val.Run.Var(Field.t.ofInt32 i)
          incr pos
        let readQ () = failwith "CLI field inputs not suported yet"
        VM.Code.iter readI readQ t 
        Field.IO.fclose h
        command args

  | "--input"::name::args ->
        // set the contents of the bank to command-line integers
        let b = bank name 
        let size = b.spec.Length 
        let contents, args = parseArgs size args
        for i = 0 to size-1 do
          b.value.[i] <- Val.Run.Var(Field.t.ofInt32 contents.[i]) 
        command args

  | "--output"::name::"to"::filename::args ->
        let b = bank name
        let t = LLVM.Tvar("%struct." + "bank_" + name)
        let h = Field.IO.fopen filename Field.IO.Write
        let format (v: Val.Run.v) = 
          match v.Arith with 
          | Val.Run.Var a ->
              // saving as a signed int, which may require renormalisation 
              try int a.Big with _ -> int (a.Big - Field.two32)
          | _ -> failwithf "format %A" v
        let pos = ref 0
        let writeI size =
          assert(size=32) 
          let i = format b.value.[!pos]
          Field.IO.write_int h i
          incr pos
        let writeQ () = failwith "non suported yet"
        VM.Code.iter writeI writeQ t 
        Field.IO.fclose h
        command args

  | "--output"::name::args -> // print the current contents of the bank 
        let b = bank name
        printfn "bank %s = %A" b.name b.value 
        command args

  | "--share"::name::args ->
        let b = bank name
        b.shared <- true
        command args

  | "--schedule"::file::args ->
        Bank.TraceReader.read file
        command args

  | "--compile"::args when Lib.w -> command (Map.fold (fun args _ name -> compile args name) args Bank.TraceReader.functions)
  | "--compile"::name::args      -> command (compile args name)

  | "--partition"::graphfile::np::name::args -> // automated QAP partitioning (experimental)
        let _, formals, body = VM.Code.findFunDef ("@"+name)
        let actuals,cfgs,args = 
          List.fold 
            (fun (actuals,cfgs,args) (t,name) -> 
              match t with 
              | LLVM.Ptr(LLVM.Tvar n) when n.StartsWith "%struct.bank_" 
                    -> let b = bank (n.Substring ("%struct.bank_".Length)) 
                       actuals @ [name,Bank.Bank b], cfgs, args
              | LLVM.Ptr(t)
                    -> // we need concrete values for t, for now only on the command line 
                       let size = VM.Code.tSize t
                       let vals, args = parseArgs size args
                       actuals @ [name,Bank.Cfg(vals)], cfgs @ [vals], args
              | _ -> failwithf "argument type %A" t
            )
            ([],[],args)
            formals
        let f = List.find (fun f -> f.name = name && f.args = cfgs) s.functions
        
        QAP.load f.qap
        let parts = System.Int32.Parse np
        Lib.time "graph extraction" (QAP.metis graphfile) parts
        command args 
        
  | "--bootstrap"::args -> Bank.bootstrap <- true; command args

  | "--keygen"::args ->
        QAP.clear() 

        // for simplicity, we take as argument the list of banks shared between proofs (we might infer it)
        // TODO: avoid extra variables when the bank is "shared" in a single mQAP 
        let shared_names, args = parseStrings args
        Array.iter (fun n -> s.banks.[n].shared <- true) shared_names  
        let summary = Lib.time "Key generation" keygen ()
        printfn "Multi-QAP with %s" summary

        match Lib.save with
        | Some mode ->
            let name = if Bank.bootstrap then "cVerifKey_BOOT.0" else "cVerifKey.0" 
            Proof.Verify.saveKey mode name Bank.s.vk.Value  // CLI bootstrapping
        | None -> ()

        command args

  | "--prove"::args when Lib.w ->
        // no need to initialize the abstract bank state here, 
        // as the C program should include a call to init. 
        let result = Lib.time "Evaluate & prove (full run)" QAPEval.main QAPEval.Prove
        if Lib.v then printfn "Evaluation returns %d" result
        command args

  | "--verify"::args when Lib.w ->
        // no need to initialize the abstract bank state here, 
        // as the C program should include a call to init. 
        let result = Lib.time "Evaluate & verify (full run)" QAPEval.main QAPEval.Verify
        if Lib.v then printfn "Evaluation returns %d" result
        command args

  | "--eval"::name::args -> // QAP-assigning evaluation & proof generation
        let _,formals,body = VM.Code.findFunDef ("@"+name)
        let actuals,cfgs,args = 
          List.fold 
            (fun (actuals,cfgs,args) (t,name) -> 
              match t with 
              | LLVM.Ptr(LLVM.Tvar n) when n.StartsWith "%struct.bank_" 
                    -> let b = bank (n.Substring ("%struct.bank_".Length)) 
                       actuals @ [name,Bank.Bank b], cfgs, args
              | LLVM.Ptr(t)
                    -> // we need concrete values for t, for now only on the command line 
                       let size = VM.Code.tSize t
                       let vals, args = parseArgs size args
                       actuals @ [name,Bank.Cfg(vals)], cfgs @ [vals], args
              | _ -> failwithf "argument type %A" t
            )
            ([],[],args)
            formals
        let f = List.find (fun f -> f.name = name && f.args = cfgs) s.functions
        QAP.load f.qap
          
        if Lib.v then
          printfn "Evaluating %s with arguments %A" name actuals 
          printfn "%d equations on %d variables" QAP.nRoot QAP.nWire 

        // execute the program on actual inputs, with exactly the same control flow as for keygen 
        // to produce the actual outputs and intermediate wire values
        // (We have a separate counter Run.nWire, so no need to reset) 
        // We get updated named banks + local assignments                
        Lib.time "QAP evaluation" (QAPEval.eval ("@"+name) actuals) [||]
        
        // process the results, preparing the commitments
        let internals = Array.init Val.Run.inter.Count (fun i -> Val.Run.inter.[i+1])

        if Lib.h then // display statistics on the range of internal values (informing exponentiation tables)
          let histogram  = Array.zeroCreate 255; 
          Array.iter (fun (x: Field.t) -> let i = QAP.log2 (x.Big + QAP.one) in histogram.[i] <- histogram.[i] + 1) internals
          printfn "#bit #variable\n"
          Array.iteri (fun i n -> if n <> 0 then printfn "%3d %7d" i n) histogram
        
        let getvals x = Array.map (fun w -> match w with Val.Run.Var v -> v | _ -> failwith "bad %A" w) x
        let finals = Array.create f.banks.Length [||] 
        begin
          let i = ref 0
          List.iter 
            (fun (name,v) -> 
              match v with 
              | Bank.Bank b -> finals.[!i] <- getvals b.value; incr i  
              | _           -> () // compile-time
            )
            actuals
        end
        if Lib.check then // check validity of all generated equations (for sanity)
          let ws = ([| Field.t.one |] :: Array.toList finals @ [internals]) 
          if Lib.v then printfn "ws %A" ws     
          Lib.time "Checking equations"   QAP.check (Array.concat ws) 

        enc_count_reset()
        // by default, the prover commits to all shared banks;
        // an optional list of banks indicates those to be recomputed by the verifier 
        let recomputed, args = parseStrings args
        let cs = Lib.time "Generating commitments"  (Bank.commit f recomputed finals) internals
        let pi = Lib.time "Generating proof"        (Proof.Prove.proof s.ek.Value) cs
        if Lib.v then 
          printfn "Non zero encoding: %d of %d" (enc_count_nonzero_get()) (enc_count_total_get()) 

        // ad hoc CLI bootstrapping experiment 
        if Lib.save = Some Lib.Simple then 
          if name = "outsource" then
            printfn "saving files for bootstrapping two-matrix"
            // This is temporary, towards bootstrapping.
            cs.[3].Save Lib.Simple "cCommitment.3" 
            pi.Save     Lib.Simple "cProof.0" 

          if name = "inner" then 
            // printfn "saving files for bootstrapping loop"
            let np = sprintf "cProof.%d" Lib.iloop
            printfn "saving proof to %s." np
            pi.Save Lib.Simple np 
            if cs.Length = 6 
              then // shared-bank variant
                if Lib.iloop % 2 = 0 then
                  let n1 = sprintf "cCommitment.%d" (2*Lib.iloop)
                  let n5 = sprintf "cCommitment.%d" (2*Lib.iloop+1)
                  printfn "saving LOCAL_0 to %s\nsaving ODD to %s." n1 n5
                  cs.[1].Save Lib.Simple n1
                  cs.[5].Save Lib.Simple n5
                else
                  let n3 = sprintf "cCommitment.%d" (2*Lib.iloop)
                  let n4 = sprintf "cCommitment.%d" (2*Lib.iloop+1)
                  printfn "saving LOCAL_1 to %s\nsaving EVEN to %s" n3 n4
                  cs.[3].Save Lib.Simple n3
                  cs.[4].Save Lib.Simple n4
              else // public-bank variant
                let n = sprintf "cCommitment.%d" Lib.iloop
                printfn "saving LOCAL to %s\n" n
                cs.[4].Save Lib.Simple n
            Lib.iloop <- Lib.iloop + 1 
                  
        // for convenience, we verify immediately after proving.
        let publicLength = Array.fold (fun a (b:b) -> if b.shared then a else a+1) 0 f.banks
        let verifierFinals = Array.sub finals 0 publicLength
        
        let csv = Lib.time "Verifying commitments"  (commitVerify f recomputed cs) verifierFinals
        let verified = Lib.time "Verifying proof"   (Proof.Verify.proof s.vk.Value pi) csv
        if not verified then 
          printfn "  cs  = %A" cs
          printfn "  csv = %A" csv
          failwithf "Proof verification failed for %s" name
          
        command args 

  | args ->
        printfn "command-line syntax"
        printfn "all commands are executed sequentially, using current values for banks, commitments, and proofs"
        printfn "integers are either comma-separated constants or random"
        printfn "  -r                       print machine-readable timings "
        printfn "  -t                       trace mode"
        printfn "  -v                       verbose mode"
        printfn "  -vv                      extra verbose mode (all add_poly/addc)"
        printfn "  -e                       print generated QAPs (before overlapping)"
        printfn "  -h                       print #bits/variable statistics"
        printfn "  -q                       quiet mode"
        printfn "  --gc                     forces GC to get stable perfs"
        printfn "  --check                  check-all-equations mode"
        printfn "  --fflib GB n             bound FFLIB table size"
        printfn "  --fflib v n              set FFLIB verbosity"
        printfn "  --fflib fpc n            set FFLIB free precomputation"
        printfn "  --encoding ENCODING      set field encoding"
        printfn "  --format [Simple|Rich]   set file binary format (Simple for bootstrapping)"
        printfn "%s\n" Field.encoding_help
        printfn "command-line mode"
        printfn "  --load PROGRAM.s              load assembly code"
        printfn "  --compile FUNCTION integers   compile named function with those compile-time values"
        printfn "  --keygen BANKS                generate mQAP and keys for all functions compiled so far"
        printfn "  --input BANK integers         set named bank to those values"
        printfn "  --input BANK from FILE        set named bank to values read from FILE"
        printfn "  --output BANK [to FILE]       output values of named bank [into FILE]"
        printfn "  --eval FUNCTION integers      evaluate and prove FUNCTION with those compile-time values";
        printfn ""
        printfn "whole-program mode"
        printfn "  -w                            select whole-program mode"
        if args = ["--help"] 
        then 0 
        else printfn "unrecognized command: %A" args; exit -1 

[<EntryPoint>]
let Main args = 
  Field.setEncoding Field.ENCODING_ENIGMA_BN   
  Lib.CPP.set()
  //let l = System.Diagnostics.Debugger.Launch()
  let r = command (List.ofArray args)

  if false then   // some rough memory freeing, to track memory leaks
    clear s;
    for i in 0..System.GC.MaxGeneration do
      printfn "Generation %d collected %d times" i (System.GC.CollectionCount i)
    printfn "managed memory      :  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    s.functions <- []
    System.GC.Collect()
    printfn "after clearing funs :  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    s.banks <- Map.empty
    System.GC.Collect()
    printfn "after clearing banks:  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    clear s
    QAPGen.Data.init()  // clearing the cache
    QAPEval.Data.init() // clearing the cache
    printfn "after clearing state:  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    QAP.clear()
    printfn "after clearing QAP  :  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    VM.Code.clear()
    printfn "after clearing code :  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    System.GC.Collect()
    printfn "after GC (leaked)   :  %13d bytes" (System.GC.GetTotalMemory true) //; Lib.pause()
    
  Lib.CPP.cleanup()
  if Lib.v || Lib.t then printfn "Done."
  r

