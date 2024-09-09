# HelixCalc

This program generates the coordinates of a monoatomic helix aligned with the z-axis.

To build use
```
make
```

To run use
```
bin/HelixCalc {distance} {angle} {dihedral} {# of atoms}
```

The program will generate two files output1.xyz and output2.xyz. Output1 contains the coordinates before being aligned with the z-axis, and Output2 is after. Two have been provided in this repository for the inputs 1.88519, 1.95676, 1.57023, 100
