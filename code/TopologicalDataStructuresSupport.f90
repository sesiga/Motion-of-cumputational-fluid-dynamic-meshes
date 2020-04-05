!TopologicalDataStructuresSupport
!
!contiene los tipos para crear las listas enlazadas

module TopologicalDataStructuresSupport

   type EdgeType_elements
       integer :: ID
       type(EdgeType_elements), pointer :: next => null()
   end type EdgeType_elements
   
   type EdgeType_facets
       integer :: ID
       type(EdgeType_facets), pointer :: next => null()
   end type EdgeType_facets
   
   type EdgeType
       integer :: Nelements
       integer :: Nfacets
       integer :: nodes(2)
       type(EdgeType_elements), pointer :: elements => null()
       type(EdgeType_facets), pointer :: facets => null()
       type(EdgeType), pointer :: next => null()
   end type EdgeType
   
   type FacetType_nodes
       integer :: ID
       type(FacetType_nodes), pointer :: next => null()
   end type FacetType_nodes
   
   type FacetType_edges
       integer :: ID
       type(FacetType_edges), pointer :: next => null()
   end type FacetType_edges
   
   type FacetType_elements
       integer :: ID
       type(FacetType_elements), pointer :: next => null()
   end type FacetType_elements
   
   type FacetType
       integer :: Nelements
       integer :: Nnodes
       integer :: Nedges
       type(FacetType_nodes), pointer :: nodes => null()
       type(FacetType_edges), pointer :: edges => null()
       type(FacetType_elements), pointer :: elements => null()
       type(FacetType), pointer :: next => null()
   end type FacetType
   
   type ElementType
       logical, pointer :: edges_state(:)
   end type ElementType
   
   type NodeType_edges
       integer :: ID
       type(NodeType_edges), pointer :: next => null()
   end type NodeType_edges
   
   type NodeType_facets
       integer :: ID
       type(NodeType_facets), pointer :: next => null()
   end type NodeType_facets
   
   type nodeType_elements
       integer :: ID
       type(nodeType_elements), pointer :: next => null()
   end type nodeType_elements
   
   type NodeType
       integer :: Nedges
       integer :: Nfacets
       integer :: Nelements
       type(NodeType_edges), pointer :: edges => null()
       type(NodeType_facets), pointer :: facets => null()
       type(Nodetype_elements), pointer :: elements => null()
   end type NodeType

end module TopologicalDataStructuresSupport