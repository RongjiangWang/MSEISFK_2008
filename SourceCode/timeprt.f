	subroutine timeprt(nt)
	implicit none
c
	integer nt
c
	integer h,m,s
c
	h=nt/3600
	m=(nt-h*3600)/60
	s=nt-h*3600-m*60
	write(*,'(a$)')'The run time = '
	if(h.gt.9)then
	  write(*,'(i5,a$)')h,':'
	else if(h.gt.0)then
	  write(*,'(a,i1,a$)')'   0',h,':'
	else
	  write(*,'(a$)')'   00:'
	endif
	if(m.gt.9)then
	  write(*,'(i2,a$)')m,':'
	else if(m.gt.0)then
	  write(*,'(a,i1,a$)')'0',m,':'
	else
	  write(*,'(a$)')'00:'
	endif
	if(s.gt.9)then
	  write(*,'(i2)')s
	else if(s.gt.0)then
	  write(*,'(a,i1)')'0',s
	else
	  write(*,'(a)')'00'
	endif
c
	return
	end
