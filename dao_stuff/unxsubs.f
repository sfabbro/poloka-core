      subroutine fabort
c      i = ieee_handler ('clear', 'all', ' ')
      return
      end!
c
c=================================================================
c
      subroutine clfile (lun)
      close (lun)
      end!
c
c=========================================================================
c
      subroutine byebye
c      character*80 out
c      call ieee_flags ('clear', 'exception', 'underflow', out)
c      call ieee_flags ('clear', 'exception', 'inexact', out)
      write (6,*)
      write (6,*) 'Good bye.'
      write (6,*)
      stop
      end!
c
c==========================================================================
c
      subroutine  oops
c      character*80 out
c      call ieee_flags ('clear', 'exception', 'underflow', out)
c      call ieee_flags ('clear', 'exception', 'inexact', out)
      write (6,*)
      write (6,*) 'Sorry about that.'
      write (6,*)
      stop
      end!
c
c==============================================================================
c
      subroutine infile (lun, file, istat)
C      character*132 expand
      character*256 expand
      character*(*) file
      if (file .eq. ' ') go to 999
  100 if (file(1:1) .eq. ' ') then
         file = file(2:len(file))
         go to 100
      end if
      open (lun, file=expand(file), status='old', err=999)
      istat = 0
      return
  999 istat = -1
      return
      end!
c
c========================================================
c
      character*(*) function expand(file)
      character*(*) file
      k = 0
      if ((file(1:1) .eq. '$') .or. (file(1:1) .eq. '~')) then
         do i=2,29
            if (file(i:i) .eq. '/') then
               call getenv (file(2:i-1), expand)
               do j=1,len(expand)
                  if (expand(j:j) .ne. ' ') k=j
               end do
               if (k .ne. 0) then
                  expand = expand(1:k)//'/'//file(i+1:len(file))
                  return
               end if
            end if
         end do
      end if
c
      do i=2,29
         if (file(i:i) .eq. ':') then
            call getenv (file(1:i-1), expand)
            do j=1,len(expand)
               if (expand(j:j) .ne. ' ') k=j
            end do
            if (k .ne. 0) then
               expand = expand(1:k)//'/'//file(i+1:30)
               return
            end if
         end if
      end do
      expand = file
      return
      end!
c
c======================================================
c
      subroutine outfil (lun, file, istat)
C      character*132 expand
      character*256 expand
      character*(*) file
      character*256 cfile
C      character*30 answer
      character*256 answer

      logical exist
      if (file .eq. ' ') then
         istat = -1
         return
      end if
 1000 if (file(1:1) .eq. ' ') then
         file = file(2:len(file))
         go to 1000
      end if
 2000 inquire (file=expand(file), exist=exist)
      if (exist) then
         cfile = file
C         call stupid ('This file already exists: '//cfile)
         answer = 'OVERWRITE'
C         write(6,*) '   Overwriting '//cfile
C         call getnam ('New output file name:', answer)
         if (answer .eq. 'OVERWRITE') then
            open (lun, file=expand(file), status='old')
            close (lun, status='delete')
         else if (answer .eq. 'END-OF-FILE') then
            call tblank
            call stupid ('Appending to existing file.')
            open (lun, file=expand(file), status='old', 
     .            iostat=istat)
            file = 'APPEND'
 1100       read (lun,110,end=9000) answer
  110       format (a)
            go to 1100
         else
            file = answer
            go to 2000
         end if
      end if
      open (lun, file=expand(file), status='new', iostat=istat)
 9000 return
      end!
c
c====================================================
c
      subroutine  delfil  (lun, file, istat)
      character*256 expand
      character*(*) file
      open (lun, file=expand(file), status='old', iostat=istat)
      if (istat .eq. 0) close (lun, status='delete')
      return
      end!
c
c====================================================
c
      character*(*) function case (string)
      character*(*) string
c
c For UNIX, leave the cases of the characters alone!
c
      case = string
      return
      end!
c
c=======================================================================
c
      subroutine ovrwrt (line, iwhich)
      character*(*) line
      if (iwhich .eq. 1) then
         write (6,1) line
    1    format (/a)
      else if (iwhich .eq. 2) then
         write (6,2) line, char(13)
    2    format (a, a1, $)
      else if (iwhich .eq. 3) then
         write (6,3) line
    3    format (a)
      else
         write (6,4) line, char(13)
    4    format (/a, a1, $)
      end if
      call flush (6)
      return
      end!
