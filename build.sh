mpif90 -c DNSbc.f90
mpif90 -c DNSbcFunctions.f90
mpif90 DNSbcTest.f90 DNSbc.o DNSbcFunctions.o -o Test
ar rcs libdnsbc.a DNSbc.o DNSbcFunctions.o
