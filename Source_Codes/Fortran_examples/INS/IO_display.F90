subroutine IO_display()

         use Grid_data

         implicit none

         if (Lx == 1.0) then

                 print *,"Uniform Grid"

         else if (Lx == 2.0) then
         
                 print *,"Non Uniform Grid"

         end if


        if (Ly == 2.0) then

                  print *,"Poisson Solver is Jacobi "
        
         else if(Ly == 3.0) then

                  print *,"Poisson Solver is GS "

          else if (Ly == 4.0) then


                   print *,"Poisson Solver is GSOR"

       end if 

end subroutine IO_display
