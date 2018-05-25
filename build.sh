mpif90 -O3 -c DNSbc.f90 -mcmodel=large -Wall -Wextra
mpif90 -O3 -c DNSbcFunctions.f90 -mcmodel=large -Wall -Wextra
mpif90 -O3  DNSbcTest.f90 DNSbc.o DNSbcFunctions.o -o Test -mcmodel=large -Wall -Wextra
ar rcs libdnsbc.a DNSbc.o DNSbcFunctions.o
