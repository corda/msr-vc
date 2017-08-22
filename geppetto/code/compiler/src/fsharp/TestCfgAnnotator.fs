module TestCfgAnnotator

open CfgAnnotator

// A simple graph structure with unique integers at each node
type DebugGraph (adj : Map<int, int list>) = 
    let mutable pairs : (int*int) list = []
    member x.ifElseList = pairs
    interface int Graph with
        member g.nodes = List.map fst <| Map.toList adj
        member g.outEdges node = adj.[node]
        member g.markIfElse st en = pairs <- (st,en)::pairs

// Executes a given function and expects it to fail (i.e. throw FailureException)
let assertFail (f:unit->unit) : unit =
    try f(); assert false
    #if DEBUG
    with Failure msg -> printfn "Case passed: %s" msg
    #else
    with Failure _ -> ()
    #endif

let goodGraphBranches : (int*int) list ref = ref []

let randomGoodGraph (n:int) : DebugGraph =
    let r = new System.Random()
    goodGraphBranches := []
    let nodeCount = ref 0 
    let res = ref <| Map<int,int list> []
    let newNode() = 
        let id = !nodeCount
        incr nodeCount
        res := (!res).Add(id,[])
        id
    let addEdge st en = res := (!res).Add(st,en::(!res).[st])
    let rec randomConnect n st en canBranch = 
        if n = 0 then 
            addEdge st en
        else
            let branchi = if canBranch then 1 else 0
            match r.Next(0,2+branchi) with
            | 0 ->  let mid = newNode()
                    let n1 = r.Next(0,n)
                    randomConnect n1 st mid canBranch
                    randomConnect (n-1-n1) mid en true
            | 1 ->  let n1 = r.Next(1-branchi, n+1)
                    if n1 > 0 then
                        let mid = newNode()
                        randomConnect (n1-1) st mid canBranch
                        randomConnect (n-n1) mid st false
                        addEdge mid en
                        if st=en then goodGraphBranches := (mid,en)::!goodGraphBranches
                    else 
                        assert canBranch
                        randomConnect n st st false
                        addEdge st en
            | _ ->  let n1 = r.Next(0,n+1)
                    assert canBranch
                    randomConnect n1 st en false
                    randomConnect (n-n1) st en false
                    goodGraphBranches := (st,en)::!goodGraphBranches
                
            
    if n = 1 then ignore <| newNode()
    else randomConnect (n-2) (newNode()) (newNode()) true
    DebugGraph(!res)

let handWrittenTests() = 
    let g = DebugGraph(Map<int,int list> 
                [ 0, [1]
                  1, [2]
                  2, [3;5]
                  3, [4;5]
                  5, [1]
                  4, []
                ])
    assertFail <| fun () -> markBranches g

    let g = DebugGraph(Map<int,int list>
                [ 0, [1;2]
                  1, []
                  2, [3;4]
                  3, [4]
                  4, [0]
                ])
    markBranches <| g
    assert (g.ifElseList = [2,4])

    let g = DebugGraph(Map<int,int list> 
                [ 0, [3]
                  3, [4]
                  4, [2;3]
                  2, [5;6]
                  6, [5;7]
                  7, [5;6]
                  5, [8;1]
                  8, [2;9]
                  9, [2;8]
                  1, []
                ])
    markBranches <| g
    assert (List.sort g.ifElseList = [2,5; 6,5; 8,2])


    let g = DebugGraph(Map<int,int list> [ 1, [2]; 2, [1;3]; 3, []])
    markBranches g
    assert (g.ifElseList = [])

    let g = DebugGraph(Map<int,int list> [ 1, [4]; 4, [2;5]; 2, [1;3]; 3, [1]; 5, []])
    markBranches g
    assert (g.ifElseList = [2,1])

    let g = DebugGraph(Map<int,int list> [ 1, [2;3]; 2, [4;5]; 3, [4]; 4, [6]; 5, [6]; 6, []])
    assertFail <| fun () -> markBranches g

    let g = DebugGraph(Map<int,int list> [ 1, [2;3]; 3, [2]; 2, [4]; 4, [3;5]; 5, []])
    assertFail <| fun() -> markBranches g

    let g = DebugGraph(Map<int,int list> 
                [ 0,  [1]; 1, [2;3]
                  2,  [4]; 3, [4]
                  4,  [5;6]
                  5,  [9]; 6, [7]
                  9,  [10]
                  10, [7;9]
                  7,  [1;8]
                  8,  [11;1]
                  11, []
                ])
    markBranches g
    assert (List.sort g.ifElseList = [(1,4);(4,7)])

    let g = DebugGraph(Map<int,int list>
                [0, [1]; 1, [2]; 2, [3]
                 3, [1;4]
                 4, [5;2]
                 5, []
                ])
    assertFail <| fun () -> markBranches g


let randomTests graphCount graphSize =
    for i = 1 to graphCount do 
        let g = randomGoodGraph graphSize
        let expected = new Set<int*int> (!goodGraphBranches)
        markBranches g
        let obtained = new Set<int*int> (g.ifElseList)
        // I wish this was an equality check, but sometimes we create
        // branches without realizing it
        assert (expected.IsSubsetOf obtained)


let llvmPrintIfElse () = 
    List.iter (function
                | LLVM.FunDef ((_,name,_),blocks) -> 
                    printfn "Processing function %s" name
                    let cfg = new LLVMCfg(blocks)
                    markBranches cfg
                    List.iter (fun (st,en) -> printfn "%s %s" st en) cfg.ifElseList
                | _ -> ()
              ) VM.code