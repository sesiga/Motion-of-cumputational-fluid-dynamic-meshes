!**************************************************
!Scatter File Data Class
!**************************************************
!Description:
!This module contains different types of Scatter data and its attributes
!
!-----------------------------------------------------------------------

module ScatterFileDataClass

   !------modules used-------------

   use precision, only : wp
   
   use constants, only : length 
   
   use utilities, only : reporterror
   
   !-----end modules used------------
   
   implicit none
   
   !Local parameters
   character(*), parameter, private :: MODname = 'ScatterFileDataClass'
   
   !Global parameters
   type ScatterDataSet
	   integer :: id !id del DataSet (unit file)
	   integer :: nvars !numero de variables 
   end type ScatterDataSet
   
   type ScatterDataFile
	   integer :: no_of_points !numero de puntos de la malla
	   real(wp), pointer :: points(:,:) !puntos con deformaion conocida
	   real(wp), pointer :: values(:,:) !deformacion de los puntos
	   integer, pointer :: global_id(:) !point location in the array
   end type ScatterDataFile
   
   type ScatterFileClass
	   character(len=LENGTH) :: filename
	   type(ScatterDataSet) :: DataSet
	   type(ScatterDataFile) :: ScatterData
   end type ScatterFileClass
   

end module ScatterFileDataClass