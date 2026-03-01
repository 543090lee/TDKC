# Discrete-Distilled-Model

writing in rust

TODO: 

add an option/parameter where if you do --consensus/no-fp then when we build, we look at accession information, and based on how well it's distributed or hitting many same targets, we include it or not.
Good idea

We dont need to index every sequence in the fasta file. only the ones that have target kmers.

~~make load and writing bulk, minimize write_all : loead is now fixed, maybe faster write~~

~~i might be using 2x memory when loading Accession registry since it's making id_to_name and name_to_id. but i dont need name_to_id during query time.~~

~~paired end read feature: just concatenate and make it coverage for that one long string~~

~~it probably has to be streamlined i think, for writing output file. kinda inefficient i think~~

Think about how fulgor does accession/color tracking, since in case if a majority of kmers are conserved then, RIP

~~SIMD vectorization~~

~~stream FASTQ not load all at once~~

reasoning on why to use 30% coverage, maybe compare it to kraken2 how they look at actual minimizer sampling

mmap

critical error in database building, like k-mer level dedup which meant if the same k-mer appeared in sequences A and B, only A's accession made it. 

CSR is kinda bad : 3.6G, so use equivalence class , delta+ VByte encoding, and maybe cite Pibiri, G. E., & Venturini, R. (2019). Techniques for Inverted Index Compression.

bro i dont know what to do when a minimizer gets two different taxid when extracting from kraken2 output