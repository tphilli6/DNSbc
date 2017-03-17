mpif90 -c -g DNSbc.f90
mpif90 -c -g DNSbcFunctions.f90
mpif90 -g DNSbcTest.f90 DNSbc.o DNSbcFunctions.o -o Test
ar rcs libdnsbc.a DNSbc.o DNSbcFunctions.o
