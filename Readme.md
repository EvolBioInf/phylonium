# Phylonium - fast and accurate estimation of evolutionary distances

This is the `phylonium` program for estimating the evolutionary distances between closely related genomes. It is much faster than alignment based approaches for phylogeny reconstruction and usually more accurate than competing alignment-free methods.


# Dependencies, Installation and Usage

This program depends on an external library: [libdivsufsort](https://github.com/y-256/libdivsufsort). It should be available for installation through a package manager of your choice. Furthermore, to build from the git repository the autotools are required.

Assuming all prerequisites are installed, the build can be started as follows. See the manual in case of compilation errors or when you try to compile on Arm: [manual](documentation/manual.pdf).

    $ autoreconf -fi -Im4
    $ ./configure
    $ make
    $ make install

After a successful build the `phylonium` executable is found in the `src` directory. It can then be run as a simple command line tool. All the sequences in one FASTA file are considered to be contigs of the same genome. The filename without the extension is used as ID in the output.

    $ phylonium Seq1.fasta Seq2.fasta
    2
    Seq1     0.0  0.1
    Seq2     0.1  0.0

The output is a distance matrix in PHYLIP format. Use `phylip neighbor`, `mat nj` from [mattools](https://github.com/kloetzl/mattools) or any other neighbor-joining implementation to build the phylogenetic tree.

For a more detailed description see the [manual](documentation/manual.pdf) or [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz903/5650408?rss=1).


# Citation

If you find this software useful, please cite our paper.

> Fabian Klötzl, Bernhard Haubold; Phylonium: Fast Estimation of Evolutionary Distances from Large Samples of Similar Genomes; Bioinformatics, Volume 36, Issue 7, 1 April 2020, Pages 2040–2046, https://doi.org/10.1093/bioinformatics/btz903


# License

Copyright © 2018 - 2025 Fabian Klötzl  
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

Individual files may be licensed differently.


# Contact

In case of bugs or unexpected errors don't hesitate to send me a mail: fabian-phylonium@kloetzl.info
