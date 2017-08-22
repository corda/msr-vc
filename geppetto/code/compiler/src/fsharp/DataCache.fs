module DataCache

// generic, functional cache: init, read, write, evict(2)
// we evict entries when the cache limit is exceeded, the guard is not longer live, or the address is no longer live
// the whole purpose is to cache the multiplications needed as we conditionally read and write, could use something less specific
// two independent instances? also, could we use one cache per gid? 
//
// TODO remove unnecessary fields

type gid = int // guard identifier
type a = Val.Native.a
type time = int

let cacheLimit = 100000

type ^v entry = { value : ^v; lastAccess : time } //$ either -1, or lastAccess recorded in lastUse! 
type key = gid * a


// knownTrue values are not stored since it is a read cache
type ^v cache = { entries        : Map<key,^v entry> // (gid,a) --> read value
                  entriesForGuard: Map<gid,Set<a>>   // gid --> { a   | (gid,a) is an entry id }
                  entriesForAddr : Map<a,Set<gid>>   //   a --> { gid | (gid,a) is an entry id }    
                  lastUse        : Map<time,key>     // last time (gid,a) was accessed (used only for eviction)
                  time           : time              // incremented at each access

                  check0         : a -> bool // used to check if the content has been initialized
                  read0          : a -> ^v 
                  write0         : a -> ^v -> unit }

let inline init isInit rd wr =  { entries         = Map.empty
                                  lastUse         = Map.empty
                                  time            = 0
                                  entriesForGuard = Map.empty
                                  entriesForAddr  = Map.empty
                                  check0 = isInit
                                  read0  = rd
                                  write0 = wr }

let  get_min m = Map.pick (fun k v -> Some (k,v)) m

let findSet k m = match Map.tryFind k m with Some v -> v | None -> Set.empty
let inline entriesForGuard g c = findSet g c.entriesForGuard
let inline entriesForAddr  a c = findSet a c.entriesForAddr

// Can also be called to initialize accesstime on first use, when e.lastAccess is -1
let inline private recordAccess (id:key) (c:^v cache) : ^v cache = 
    let e = c.entries.[id]
    let t = e.lastAccess
    let t'= c.time + 1
    { c with lastUse = c.lastUse.Remove(t).Add(t',id) 
             entries = c.entries.Add(id,{e with lastAccess=t'}) // update
             time    = t' }


let inline private evictEntry ((g,a):key) (c:^v cache) : ^v cache =
    let t = c.entries.[(g,a)].lastAccess
    let eg = c.entriesForGuard
    let ea = c.entriesForAddr
    { c with entries = c.entries.Remove(g,a)
             lastUse = c.lastUse.Remove(t)
             entriesForGuard = eg.Add(g,eg.[g].Remove(a))
             entriesForAddr  = ea.Add(a,ea.[a].Remove(g)) }

let inline private evictOldest (c:^v cache) : ^v cache =
    let id = snd (get_min c.lastUse)
    evictEntry id c

/// remove all entries for that guard
let inline evictByGuard g c =
    let aset = entriesForGuard g c
    Set.fold (fun c a -> evictEntry (g,a) c) c aset

/// remove all entries for that address
let inline evictByAddr a c =
    let gset = entriesForAddr a c
    Set.fold (fun c g -> evictEntry (g,a) c) c gset

let inline private shrinkBack c = 
    let cr = ref c
    while (!cr).entries.Count > cacheLimit do
        cr := evictOldest !cr
    !cr

let inline private newEntry ((g,a) as id:key) (v:^v) (c:^v cache) : ^v cache = 
    let e = { value = v; lastAccess = -1 }
    let aset = entriesForGuard g c
    let gset = entriesForAddr a c
    let c = 
      recordAccess id
        { c with entries = c.entries.Add(id,e) // extend
                 entriesForGuard = c.entriesForGuard.Add(g,aset.Add(a))
                 entriesForAddr  = c.entriesForAddr.Add(a,gset.Add(g)) }
    shrinkBack c

let inline read (cond:^v when ^v :(member knownTrue : bool)) 
                (gid:gid) (a:a) (c:^v cache) : (^v * ^v cache)
    =
    if (^v : (member knownTrue : bool)(cond)) 
    then 
      c.read0 a, c
    else 
      match c.entries.TryFind(gid,a) with
      | Some e -> e.value, recordAccess (gid,a) c
      | None   -> let v = cond * c.read0 a in 
                  v, newEntry (gid,a) v c

let inline write (cond:^v when ^v : (member knownTrue : bool))
                 (gid:gid) (a:a) (v:^v) (c:^v cache) : ^v cache
    =
    let c = evictByAddr a c
    if (^v : (member knownTrue : bool)(cond)) 
    then 
      c.write0 a v; c
    else 
      // insert an entry in the cache
      let v' = 
        if c.check0 a then let v0 = c.read0 a in v0 + cond*(v-v0) 
        else v
      c.write0 a v'
      newEntry (gid,a) v' c