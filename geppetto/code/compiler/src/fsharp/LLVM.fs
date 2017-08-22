module LLVM

/// Simplified abstract syntax and operations for LLVM code 
/// See LLVM reference manual for their intended semantics

type fvar = string 
type var = string
type label = string

type typ = 
  | Int of int // number of bits, in { 1 8 16 32 64 }
  | Array of int * typ
  | Struct of typ list 
  | Ptr of typ
  | Tvar of var
  | TVoid

type value = 
  | Var of var
  | Const of int
  | Computed of typ * value * value list // constant computed as getelementptr on a global
  | Void
  
type rhs = 
  | Phi of List<label*value option> 
  | BitCast of typ * var
  | Select of value * value * value 
  // memory
  | Addr  of typ * value * value list // aka getelementptr
  | Alloc of typ 
  | AllocN of typ * value
  | Load of value  // var*
  | Call of fvar * value list 
  // arithmetic  
  | Add of value * value
  | Sub of value * value
  | Mul of value * value
  | Srem of value * value
  | Eq  of value * value 
  | Neq of value * value 
  // booleans
  | And of value * value
  | Or of value * value
  | Xor of value * value
  | Shl of value * value
  | Shr of value * value

  | Sge of value * value 
  | Sle of value * value
  | Sgt of value * value 
  | Slt of value * value 
  | Ult of value * value
  | Ugt of value * value
  | Zext of value 
  
type instr = 
  | Let of var * rhs        // local assignment
  | Store of value * value  // value* := value 
  // terminators:
  | Branch of label
  | BranchIf of var * label * label // conditional branching
  | Switch of var * label * (value * label) list // switching; the first label is the default (?)
  | Return of value
 
type block = { label: label;  instr: instr list } // todo (performance): use instruction arrays instead
let mkBlock label instr = { label=label; instr=instr }

type decl =
  | GlobalDef of var * typ * int option // a bit ad hoc
  | GlobalStr of var * string // a bit ad hoc, could merge with GlobalDef
  | FunDef of (typ * string * (typ * var) list) * block list
  | TypeDef of var * typ 
