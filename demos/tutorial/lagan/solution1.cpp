#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>


//![main]
int main() {
    // read sequences from fasta files

    // create file objects for the fasta files
    CharString filePath1 = getAbsolutePath("example_reference.fa");
    CharString filePath2 = getAbsolutePath("example_query.fa");
    SeqFileIn sfi1;
    if (!open(sfi1, toCString(filePath1))) {
        std::cerr << "ERROR: Could not open the first file" << std::endl;
        return 1;
    }
    SeqFileIn sfi2;
    if (!open(sfi2, toCString(filePath2))) {
        std::cerr << "ERROR: Could not open the second file" << std::endl;
        return 1;
    }

    // read the first sequences from each file object
    CharString meta;
    Dna5String seqH;
    Dna5String seqV;
    readRecord(meta, seqH, sfi1);
    readRecord(meta, seqV, sfi2);

    //![solution]
    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;

    TIndex index(seqH);
    resize(indexShape(index),qgramSize);
    //![solution]

    return 0;
}
//![main]