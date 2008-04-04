module Nodes = AllDirNode.AllDirF
module Edges = Edge.LazyEdge
module TreeOps = AllDirChar.F
module CharOps = AllDirChar.CharScripting

module M = Scripting.Make (Nodes) (Edges) (TreeOps) (CharOps)

open M
include M
