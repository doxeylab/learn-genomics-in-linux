# Task1 - Learning the Linux shell

This task will introduce students to the Linux command-line (shell environment).

## What is the Linux Shell?

The <i>shell</i> is a command-line programming language for interacting with the UNIX operating system.  https://en.wikipedia.org/wiki/Unix_shell

There are several different shell languages. What we will be using in this course is a popular shell flavor called <b>BASH</b>

https://en.wikipedia.org/wiki/Bash_(Unix_shell)

We will be learning basic commands, but BASH is actually a language that can perform complex programming tasks.


## Getting Started

Before you can begin with the coding exercises, you must have access to a linux machine.
You can either use your own local system or a remote VM that has been set up for you.

### Accessing the remote VM

If you have been given access to a remote VM, you can access it either through the url that has been given to you, or directly through the Terminal like this:

```
ssh -i /path/to/your/.ssh/publickey yourUserName@remoteIP
```

When you are done, you can leave your session by typing...

```
exit
```

## Once you are logged on | Learning Unix

If you have logged in correctly, you should see a welcome screen.

You are now in a unix command-line environment. For more information and instructions:

* [CodeAcademy](https://www.codecademy.com/learn/learn-the-command-line) - Learn the command line
* [TeachingUnix](http://www.ee.surrey.ac.uk/Teaching/Unix/) - Another Unix Tutorial
* [BasicCommands](http://mally.stanford.edu/~sr/computing/basic-unix.html) - List of common commands

## A unix primer

### Navigating files/folders

When you log in, by default you start in your home directory.

Type...

```
pwd
```

And this will print the current location (e.g., /home/username)

You can always get back to your home folder by typing

```
cd
```

To view the contents of your current folder type

```
ls
```

To make the folder "task1", type

```
mkdir task1
```

Change directory into 'task1' folder

```
cd task1
```

And now create a file called file.txt

```
>file.txt
```

Open the file with `nano`. This is one built-in text editor. There are others.

```
nano file.txt
```

Now enter a few lines of text, type 'ctrl-o' and then 'ctrl-x' to save and exit

To print the contents of your file type

```
cat file.txt
```

To count the number of words and lines in your file

```
wc file.txt #words
wc -l file.txt #lines
```

To copy the file to a new file

```
cp file.txt newfile.txt
```

You can also 'move' a file to a new location or rename it using `mv`.

To combine both files together into a third file

```
cat file.txt newfile.txt > thirdfile.txt
```

The '>' redirects the output of the commands on the left of it to a file specified on the right.

Delete the file

```
rm file.txt
```

Move back to the previous folder

```
cd ..
```

And delete the folder 'task1'

```
rmdir task1
```

### Getting help on linux commands and program usage

For most commands, you can get more information on their usage by typing `man` 'command'.

e.g., try

```
man ls
#type 'q' to quit
```

`man` will work with some bioinformatics tools. However, not always.

e.g., for help on the `blastp` tool, type

```
blastp -h # or
blastp --help
```

### Additional operations

#### Pattern finding with grep

```
grep "word" file.txt  # prints lines in file.txt containing "word"
```

We can also chain together multiple commands like this using the `|` (pipe) operator.

```
grep "word" file.txt | wc -l  # will count the number of lines containing the word "word"

# or alternatively

cat file.txt | grep "word" | wc -l  # does the same thing as above
```

#### copying a file to and from a remote server

To
```
scp /path/to/file.txt username@remoteserver.com:/path/to/location/.
```

From
```
scp username@remoteserver.com:/path/to/file.txt /path/to/location/.
```

#### downloading a file from the internet

```
wget <url>
```

#### file compression/uncompression

This is done using programs such as `tar`, and `gzip` and `gunzip`.
e.g.,
```
gzip file.txt # to compress it
gunzip file.txt.gz # to uncompress it
```


# ASSIGNMENT QUESTIONS

Your assignment will be to write a series of shell commands to answer the following questions. Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.

Hint: remember to use `man` if you want to explore added functionality of commands.

* Download this file and extract its contents https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/e-coli-k12-genome.fasta.gz

* How many characters are in this file?

* How large is the file in megabytes?

* 


* Download


