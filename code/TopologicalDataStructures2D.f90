!TopologicalDataStructures2D
!
!contiene los tipos utilizados en la descripción de 
!la malla en 2D
    
module TopologicalDataStructure2D

   !modules used
   use precision, only : wp
   
   use constants, only : LENGTH
   !end modules used

   type ElementType
       character(LENGTH) :: figure
       integer :: ID
       integer :: Nnodes !numero de nodos
       integer :: Nedges !numero de aristas
       
       integer, pointer :: nodes(:) !nodos del elemento
       integer, pointer :: edges(:) !aristas del elemento
       type(ElementType), pointer :: next => null()
   end type ElementType
   
   type NodeType
       integer :: ID
       integer :: Nelements !numero de elementos
       integer :: Nedges !numero de aristas
       
       integer, pointer :: elements(:) !elementos que contienen ese node
       integer, pointer :: edges(:) !aristas que contienen ese node
       type(NodeType), pointer :: next => null()
   end type NodeType
   
   type EdgeType
       integer :: ID
       integer :: Nelements !numero de elementos
       integer :: Nnodes !numero de nodos
       
       integer, pointer :: elements(:) !elementos que contienen esa arista
       integer :: nodes(2) !nodos que contienen esa arista
       type(EdgeType), pointer :: next => null()
   end type EdgeType
   
   type TopologyType
       integer :: Nnodes !numero total de nodos
       integer :: Nelements !numero total de elementos
       integer :: Nedges !numero total de aristas
       
       type(ElementType), pointer :: elements(:)
       type(NodeType), pointer :: nodes(:)
       type(EdgeType), pointer :: edges(:)
   end type TopologyType
   
   type NodeTypeElements
       integer :: element
       type(NodeTypeElements), pointer :: next => null()
   end type NodeTypeElements
   
   type EdgeTypeSupport
       integer :: element
       type(EdgeTypeSupport), pointer :: next => null()
   end type EdgeTypeSupport
   
end module TopologicalDataStructure2D