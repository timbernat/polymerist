'''Tools for manipulating, processing, and pretty-printing text from files and string-like objects'''

from string import ascii_letters, ascii_lowercase, ascii_uppercase
ascii_printable = ''.join(chr(i) for i in range(33, 127)) # ASCII printable characters, minus SPACE (" ", 32) and DELETE (127)