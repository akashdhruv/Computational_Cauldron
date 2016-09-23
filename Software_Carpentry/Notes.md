
- To download from a web page type

  ~~~terminal
  wget link_to_download_location
  ~~~

- To unzip

  ~~~terminal
  unzip file_name.zip (also look up 'tar')
  ~~~

- Print username and current directory

  ~~~terminal
  whoami
  pwd
  ~~~

- ls options

  ~~~terminal
  ls
  ls -F
  ls -a
  ls -l
  ls -h
  ls /
  ls -n
  ls --help
  ~~~

- rm options

  recursive remove

  ~~~terminal
  rm -r directory_name
  rm file_name
  ~~~
  
  Interactive remove: asks permission before deleting

  ~~~terminal
  rm -ri directory_name
  ~~~

  Deleting empty directories

  ~~~terminal
  rmdir name_of_empty_directory
  ~~~

- touch options

  ~~~terminal
  touch file_name
  ~~~

  Changes access time of the existing file. If the file does not exist, it creates an empty file

  Reasons to use touch:
  - Avoiding text editor
  - Updates time stamp

- word count : lines, words, characters

  ~~~terminal
  wc *.extension
  wc file_name
  ~~~ 

  Flags: -l (lines) -w (words) -c (characters)

- store output to another file

  ~~~terminal
  wc -l file_name > output_filename.txt
  ~~~

- head vs cat vs tail

  'head' prints only first 10 lines, 'cat' print everything
  'tail' prints the last 10 lines

   also look - less, sort, tee

- pipes
  
  ~~~terminal
  wc -l *.pdb | sort -n
  ~~~
  sort the output of first command numerically

  ~~~terminal
  wc -l *.pdb | sort -n | head -n 1
  ~~~
  sort the output of the first command and print the first line
   
  ~~~terminal
  screen
  mpirun -n 4 ./exec_file | tee output.txt
  ~~~

  add output to file

  ~~~terminal
  echo hello > hello.txt
  ~~~

  append to file

  ~~~terminal
  echo hello >> hello.txt
  ~~~

- loops

  ~~~terminal
  for filename in *.dat
  > do
  >      echo $filename
  >      echo 
  >      head -n 100 $filename | tail -n 5
  > done
  ~~~

- git status outside working directory

  ~~~terminal
  git --git-dir=/home/akash/Github/ParaSolve_AMR/.git --work-tree=/home/akash/Github/ParaSolve_AMR status
  ~~~



