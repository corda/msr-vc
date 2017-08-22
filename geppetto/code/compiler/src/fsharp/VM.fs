module VM

(* part of the interpreter that does not depend on the domain of values *)  

// we could turn loaded code into module global state (to avoid closure allocation)

open LLVM

/// global state: some parsed LLVM assembly file 
let mutable code: LLVM.decl list = [] 

module Code =
  // For now, code segments are the parsed lists and maps; consider using arrays instead
  type blocks = block list

  /// abstract, stateful code pointer within a function
  type xp(segment:blocks) = 
    let segment  = segment       // all code blocks for the function 
    let mutable label : string     = segment.Head.label // label of the current block
    let mutable block : instr list = segment.Head.instr // remaining instructions in the current block
    member xp.Segment = segment
    member xp.Label = label
    member xp.Next = 
      match block with 
      | x::xs -> block <- xs; x
      | []    -> failwith "bad block"
    member xp.Goto lbl =
      block <- (List.find (function (b:block) -> b.label = lbl) segment).instr
      label <- lbl

  open System.Collections
  let mutable private tdef : Hashtable = null // string -> typ
  let mutable private fdef : Hashtable = null // string -> (typ * var) list * block list
  let mutable private tsize: Hashtable = Hashtable() // typ -> int 
   
  let clear() = 
    tdef <- null
    fdef <- null
    tsize<- null

  // profiler says these are worth caching --- we assume ds is fixed.
  
  let findTypeDef name: typ = 
    if tdef = null then  
      tdef <- Hashtable() 
      List.iter (function TypeDef (n,t) -> tdef.[n] <- t  | _ -> ()) code
    tdef.[name] :?> typ

  let findFunDef name: typ * (typ * var) list * block list = 
    if fdef = null then  
      fdef <- Hashtable() 
      List.iter (function FunDef ((result,n,formals),blocks) -> fdef.[n] <- (result, formals, blocks) | _ -> ()) code
    if fdef.ContainsKey name
    then fdef.[name] :?> (typ * (typ * var) list * block list)
    else failwithf "function not found: %s" name

  let globalStr name = 
    let rec f = function
      | GlobalStr(n,s)::_ when n = name -> s.Substring(2,s.Length - 6).Replace("\\0A","\n") // ad hoc, removing 'c"' and '\00"'
      | _::ds                           -> f ds 
      | []                              -> failwithf "global string %s." name
    f code

  let globalInt name = 
    let rec f = function
      | GlobalDef(n,_,Some i)::_ when n = name -> Some i
      | _::ds -> f ds 
      | [] -> None
    f code
  // we compute the *interpreter* size, that is, 
  // we count one word for every integer, irrespective of their native size
  // (we used to take a multiplier n; now applied by caller)
  let rec tSize t = 
    match t with 
    | Int _  -> 1 
    | Ptr _  -> 1
    | Array(n',t) -> n' * tSize t 
    | Struct(ts) -> (List.fold (fun a t -> a + tSize t) 0 ts)  
    | Tvar name  -> if tsize.ContainsKey name 
                    then tsize.[name] :?> int 
                    else 
                      let s = tSize (findTypeDef name) in
                      tsize.[name] <- s
                      s
    | _ -> failwithf "not supported yet: %A" t

  // type-directed data traversal, usable both for loading and saving
  let rec iter readI readQ = function
    | Int sz          -> readI sz
    | Ptr _          -> readQ() // deprecated
    | Array(n,t)     -> for i = 0 to n-1 do iter readI readQ t
    | Struct(ts)     -> List.iter (fun t -> iter readI readQ t) ts
    | Tvar "%struct.Qapfp_t" -> readQ()
    | Tvar name      -> iter readI readQ (findTypeDef name)
    | (TVoid) as t   -> failwithf "not supported yet: %A" t

  let rec offset t ns : int = 
    let r = 
      match t, ns with 
      | Tvar name           , _     -> offset (findTypeDef name) ns
      | Struct (t::ts)      , 0::ns -> offset t ns  
      | Struct (t::ts)      , n::ns -> tSize t + offset (Struct ts) ((n-1)::ns)  
      | (Ptr t | Array(_,t)), n::ns -> 
            // a more general fix is needed! 
            (if t = Int 8 && n % 4 = 0 then n/4 else n) // needed for some MapReduce examples, due to LLVM aliasing of byte*-allocated data
            // n                                        // needed for some X509 examples
            * tSize t + offset t ns 
      | _, []                           -> 0
      | _                               -> failwithf "offset %A, %A" t ns
    r

  open Microsoft.FSharp.Text.Lexing
  let load (filename: string) = 
    let lexbuf = 
      let r = new System.IO.StreamReader(filename) 
      LexBuffer<_>.FromTextReader r 
    try
      code <- LLVM_grammar.start Lexer.token lexbuf
    with e -> 
      let pos = lexbuf.EndPos
      let lastToken = new System.String(lexbuf.Lexeme)
      printfn "LLVM Parsing of %s failed line %d, column %d on token: %s" filename pos.Line pos.Column lastToken 
      exit -1