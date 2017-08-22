module CfgAnnotator

// ------------ First part works on general graphs: doesn't depend on LLVM  ---------------------

type 'Node Graph when 'Node : comparison = 
    abstract nodes : 'Node list
    abstract outEdges : 'Node -> 'Node list
    abstract markIfElse : 'Node -> 'Node -> unit

type private 'Node VisitRes = HitReturn
                            // HitInVisit (loopEntry, (exitContinuation) option)
                            | HitInVisit of 'Node * 'Node VisitRes option
                            // HitVisited (exitNode)
                            | HitVisited of 'Node

type private 'Node NestingCheck = BranchEnd of int * 'Node
                                | LoopStart of 'Node
                                | LoopEnd

// Returns the transpose of a graph, i.e. the inward edges incident on a node
let mkIndegs (g:'Node Graph) : Map<'Node,'Node list> =
    let nodes = g.nodes 
    let mutable result = new Map<'Node,'Node list> (List.map (fun n -> (n,[])) nodes)
    for i in nodes do for j in g.outEdges i do result <- result.Add(j,i::result.[j])
    result

(*
Performs a depth-first-search on the control flow graph.
Invokes cfg.markIfElse on the node pairs denoting the entry
and exit nodes of an if-else structure. Invokes failwith if
the input graph could not be obtained by a simple composition 
of if/else and loop structures (i.e. we do not support break/
continue/early returns etc.)

A little caveat: any node without outgoing edges are considered
a "return". Be a little careful if your input code might contain
exceptions and what not.
 
validCode := {<empty>} | <straight-line-code>;
           | validCode validCode
           | if(...) validCode [else validCode]
           | while(...) validCode

Should run in O(max(n,k^2) log n) time, where k is the max loop 
nesting level and log n comes from various internal Set/Map lookups 
The quadratic part can probably be improved, but I don't expect k 
to be high anyway
*)

let markBranches (cfg:'Node Graph) = 
    let visited = ref <| new Set<'Node> [||]
    let invisit = ref <| new Set<'Node> [||]

    let multipleExits() = failwith "Loop has multiple exits. Currently unsupported."
    let nestingError()  = failwith "Control structures are not properly nested"

    let indegs = Map.map (fun _ l -> ref <| List.length l) <| mkIndegs cfg
    let discov  = ref <| new Map<'Node,int> [||]  // keys(discov) = visited Union invisit
    let nodeCount = ref 0
    let validPair st en child = let den = (!discov).[en]
                                den<=(!discov).[child] && den>(!discov).[st]

    // These three functions could be faster if we used append-friendly lists
    let rec discoveryRange x = 
        match x with
        | HitReturn | HitVisited _ -> Choice1Of2 -1
        | HitInVisit(root,None) -> let d = (!discov).[root] in Choice2Of2 (d,d)
        | HitInVisit(root,Some y) -> 
            let d = (!discov).[root]
            match discoveryRange y with
            | Choice1Of2(hi)    -> assert (d>=hi); Choice1Of2(d)
            | Choice2Of2(lo,hi) -> assert (d>=hi); Choice2Of2(lo,d)
    let canAppend res1 res2 = 
        match res1 with 
        | HitReturn | HitVisited _ -> false
        | _ -> match discoveryRange res1 with
               | Choice1Of2 _ -> false
               | Choice2Of2 (lo,hi) -> match discoveryRange res2 with 
                                       | Choice1Of2(hi2) | Choice2Of2(_,hi2) -> hi2<=lo
    let rec appendInVisit res1 res2 = match res1 with HitReturn | HitVisited _  -> assert false; failwith ""
                                                    | HitInVisit(x,None) -> HitInVisit(x,Some res2)
                                                    | HitInVisit(x,Some res1') -> 
                                                        HitInVisit(x,Some(appendInVisit res1' res2))

    let rec visit (node:'Node) : ('Node VisitRes * 'Node NestingCheck list) = 
        let inC = indegs.[node]
        let popBranchEnd res joinpoint branchnest nest =
            match nest, branchnest with
            | BranchEnd (x,joinpoint') :: t1, [BranchEnd (y,joinpoint'')] 
                when joinpoint = joinpoint' && joinpoint = joinpoint'' -> 
                    let t1 = if x+y = 0 then t1 else BranchEnd(x+y,joinpoint')::t1
                    cfg.markIfElse node joinpoint; res,t1
            | _ -> nestingError()

        // Checks if the recursion has finished visiting a loop (and its descendants), and is about to return from it
        let rec fixLoopHead (result,nest) = 
            match result,nest with
            | HitInVisit(x,Some inres), (LoopStart node' :: t1 | LoopEnd :: LoopStart node' :: t1) when x=node-> 
                assert (node'=node) // Can malformed graphs violate this? If so this should be nestingError()
                fixLoopHead(inres,t1)
            | HitInVisit(x,None),_ when x=node && (cfg.outEdges x).Length = 1 -> 
                failwith "Can't reach end of function; possibly reached infinite loop"
            | HitInVisit(x,_),_ when (!visited).Contains(x) -> nestingError() // invisit is suddenly visited
            | _ -> result,nest
                
        // Eventually this should go down to zero
        decr inC
        if (!visited).Contains node then HitVisited node,[BranchEnd (-1,node)]
        else if (!invisit).Contains node then HitInVisit (node,None),[LoopStart node]
        else 
            let children = cfg.outEdges node
            invisit := (!invisit).Add(node)
            discov := (!discov).Add(node,!nodeCount)
            let mediscov = !nodeCount
            incr nodeCount;
            let result, nestingChecks =
                match children with 
                | [] -> HitReturn, []
                | [x] -> visit x
                | [x;y] -> 
                    let xr,xnest = fixLoopHead <| visit x
                    let yr,ynest = fixLoopHead <| visit y
                    match xr,yr with
                    // Corner case: nested loops transformed into an if-else condition
                    | HitInVisit(xroot,None), HitInVisit(yroot,None) when xroot=yroot ->
                        assert (xnest = [LoopStart xroot])
                        assert (ynest = [LoopStart yroot])
                        cfg.markIfElse node xroot
                        xr,xnest

                    // Loop exit processing
                    | HitInVisit (xroot,_), HitVisited yend when validPair node yend y || yend = xroot -> // XXX check this again
                        popBranchEnd xr yend ynest xnest   // if-in-loop case
                    | _,_ when canAppend xr yr -> appendInVisit xr yr, [LoopEnd] @ xnest @ ynest
                    | _,_ when canAppend yr xr -> appendInVisit yr xr, [LoopEnd] @ ynest @ xnest
                    | HitInVisit _, HitInVisit _ -> nestingError()
                    | HitInVisit(_,Some _), _ 
                    | _, HitInVisit(_,Some _) -> multipleExits()

                    // If assertions fail, canAppend may be buggy
                    | HitInVisit(root,None), _ -> assert (xnest <> [LoopStart root]); nestingError()
                    | _, HitInVisit(root,None) -> assert (ynest <> [LoopStart root]); nestingError()

                    | (HitReturn | HitVisited _), HitVisited yend -> 
                        if validPair node yend y || xr = yr
                            then popBranchEnd xr yend ynest xnest // Top-level if or if-in-if case
                            else failwith "Branch never converges back"
                    | (HitVisited _ | HitReturn), HitReturn -> failwith "Early function return not supported yet"
                | _ -> failwith "Basic block has more than two children. We don't support switch-case yet."

            invisit := (!invisit).Remove node
            visited := (!visited).Add node
            let result, nestingChecks = fixLoopHead (result, nestingChecks)
            // There are three kinds of incoming edges
            let nestingChecks = if !inC = 0 then nestingChecks
                                else BranchEnd (!inC,node) :: nestingChecks
            result, nestingChecks

    in
    let assertVis (vr,_) = match vr with HitVisited _ -> ()
                                        | _ -> failwith "Entire function needs to be reachable from entry point"
#if DEBUG2
    try
#endif
    match cfg.nodes with [] -> ()
                       | x::xs -> incr indegs.[x]
                                  match visit x with HitReturn,[] -> List.iter (assertVis << visit) xs
                                                   | _ -> failwith "Function does not lead to return"
#if DEBUG2
    with ex ->  printfn "%A" !visited
                printfn "%A" !invisit
                raise ex
#endif


// ------------------------ Implements interface for LLVM -----------------------------------------------------------

open LLVM

let rec last l = match l with [] -> failwith "last of empty list"
                            | [x] -> x
                            | _::xs -> last xs

let blockEntry (b : block) = b.label
let blockExit  (b : block) = match last b.instr with
                             | Branch l -> [l]
                             | BranchIf (_,tl,fl) -> [tl;fl]
                             | _ -> []

let rec private mkMap (m : Map<string, string list>) bs = 
    match bs with
    | [] -> m
    | b::bs -> let m' = m.Add(blockEntry b, blockExit b) in mkMap m' bs

type LLVMCfg (b : block list) = 
    let adj = mkMap (Map<string,string list> [||]) b
    let mutable pairs : (string*string) list = []
    member g.ifElseList = pairs
    interface string Graph with
        member g.nodes = List.map (fun x -> x.label) b
        member g.outEdges node = adj.[node]
        member g.markIfElse st en = pairs <- (st,en)::pairs
        
    // Debugging aid
    member g.printAdj() = 
        List.iter (fun b -> printfn "%s %A" b.label adj.[b.label]) b

