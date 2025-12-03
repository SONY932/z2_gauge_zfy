Module Files_mod
   contains

     Character (len=64) function File_i( file, I)
        character (len=64) :: file
        integer            :: i
        write(File_i,'(A,"_",I0)') trim(file),i
      end function File_i

end Module Files_mod
