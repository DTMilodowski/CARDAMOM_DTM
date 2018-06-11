
subroutine oksofar(message)

! subroutine prints information to screen 
! could be modifed to ensure information 
! is passed to specific declared error files

implicit none
! local declaration
character(200), intent(in) :: message

! print the message to screen / output files
print*, message

end subroutine oksofar
