//
// Created by abadd0n on 4/30/16.
//

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;


template <typename TSeedString>
bool printSeedString(TSeedString & seedSet) {
    typedef typename Iterator<TSeedString>::Type            TSeedChainIter;
    std::cout << "SeedSet" << std::endl;
    for (TSeedChainIter it = begin(seedSet); it != end(seedSet); ++it) {
        if ( endPositionH(*it) <= beginPositionH(*it) || endPositionV(*it) <= beginPositionV(*it) )
            return false;
        std::cout << "SeqH\tStart:\t" << beginPositionH(*it) << "\tEnd:\t" << endPositionH(*it) << "\tSeqV\tStart:\t" <<
                beginPositionV(*it) << "\tEnd:\t" << endPositionV(*it) <<
                "\tupper diag:\t" << upperDiagonal(*it)<< "\tlower diag:\t" << lowerDiagonal(*it) << std::endl;
    }
    return true;
}


template <typename TSeedChain>
void updateSeedPositions(TSeedChain & seedChain, unsigned globalPosH, unsigned globalPosV)
{
    typedef typename Iterator<TSeedChain>::Type            TSeedChainIter;
    for (TSeedChainIter it = begin(seedChain); it != end(seedChain); ++it) {
        setBeginPositionH(*it, beginPositionH(*it) + globalPosH);
        setEndPositionH(*it, endPositionH(*it) + globalPosH);
        setBeginPositionV(*it, beginPositionV(*it) + globalPosV);
        setEndPositionV(*it, endPositionV(*it) + globalPosV);

        //diag
        setUpperDiagonal(*it, upperDiagonal(*it) + globalPosH-globalPosV);
        setLowerDiagonal(*it, lowerDiagonal(*it) + globalPosH-globalPosV);
    }
}

template <typename TSeedChain, typename TSeqString, typename TQGramSize, typename TSeed>
inline bool createSeedChain(TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, TQGramSize const & qgramSize, TSeed const & /*tag*/)
{
    if ( length(seqH) < qgramSize || length(seqV) < qgramSize ) {
        return false;
    }

    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;
    typedef Infix<Dna5String>::Type TInfix;

    std::clock_t start;
    start = std::clock();

    TIndex index(seqH);
    resize(indexShape(index),qgramSize);
    TInfix kmer;

    double durationqgram = (std::clock() - start )/ (double) CLOCKS_PER_SEC;
    std::cout << "Runtime build qgram : " << durationqgram << std::endl;

    //create seedSet
    typedef SeedSet<TSeed>      TSeedSet;


    TSeedSet seedSet;
    Score<int, Simple> scoringScheme( 2,-1,-2 ); // match, mismatch, gap
    //std::cout << "qPos: ";
    for (unsigned qPos = 0; qPos < length(seqV)-qgramSize+1; ++qPos)
    {
        //std::cout << qPos << std::endl;
        kmer = infix(seqV, qPos, qPos+qgramSize);
        hash(indexShape(index), begin(kmer));
        for ( auto occurence : getOccurrences( index, indexShape( index ) ) ) {
            // add seed to seed set using CHAOS chaining method
            if ( !addSeed( seedSet, TSeed(occurence, qPos, qgramSize), 5 /*max diag dist*/, 10 /*band width*/, scoringScheme, seqH, seqV, Chaos() ) )
                addSeed( seedSet, TSeed(occurence, qPos, qgramSize ), Single()); // just add seed if CHAOS fails

            // Probably we don't need to extend them manually
            //extendSeed( seed, d5s1, d5s2, EXTEND_BOTH, scoringScheme, 3 /*score drop off*/, UngappedXDrop() );
            //appendValue(seedChain, TSeed(getOccurrences(index, indexShape(index))[i], qPos, qgramSize));
        }
    }

    //std::cout << std::endl;

    if ( empty (seedSet) )
        return false;

    chainSeedsGlobally(seedChain, seedSet, SparseChaining());


    return true;
}


struct LaganOptions
{
    CharString path2file1;
    CharString path2file2;
    std::vector<unsigned> qgramsize;

    //LaganOptions() :
    //{}
};



ArgumentParser::ParseResult
parseCommandLine(LaganOptions & options, int argc, char ** argv)
{
    // Argument parser
    ArgumentParser parser("LAGAN");

    // Set short description, version, and date.
    setShortDescription(parser, "LAGAN Algorithm");
    setVersion(parser, "1.0");
    setDate(parser, "May 2016");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser,
                   "This program uses LAGAN algorithm to aligns two genomes");

    // We require one argument (fa or fq file to read sequences from).
    addArgument( parser, ArgParseArgument( ArgParseArgument::STRING, "PATH" ) );

    // User can also provide second fa or fq file in case he has genomes in separate files
    // second file is arg
    addArgument( parser, ArgParseArgument( ArgParseArgument::STRING, "PATH" ) );

    addOption( parser, ArgParseOption(
            "q", "qgramsize", "one or more qGram sizes",
            ArgParseArgument::STRING, "SIZE", true ) );

    // Parse command line.
    ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "qgramsize"))
    {
        for (std::string sizeString : getOptionValues(parser, "qgramsize"))
        {
            unsigned qsize;
            if (!lexicalCast(qsize, sizeString))
            {
                std::cerr << "ERROR: Invalid qGramSize " << sizeString << "\n";
                return seqan::ArgumentParser::PARSE_ERROR;
            } else
            {
                options.qgramsize.push_back(qsize);
            }
        }
    }

    getArgumentValue(options.path2file1, parser, 0);
    getArgumentValue(options.path2file2, parser, 1);

    return ArgumentParser::PARSE_OK;
}

int main(int argc, char ** argv) {

    // Parse the command line.
    LaganOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // read sequences from fasta files

    // create file objects for the fasta files

    CharString filePath1 = options.path2file1;
    CharString filePath2 = options.path2file2;
    if (options.path2file1[0] != '/') {
        filePath1 = getAbsolutePath(toCString(filePath1));
    }
    if (options.path2file2[0] != '/') {
        filePath2 = getAbsolutePath(toCString(filePath2));
    }
    SeqFileIn sfi1;
    if (!open(sfi1, toCString( filePath1 ) ) )
    {
        std::cerr << "ERROR: Could not open the first file" << std::endl;
        return 1;
    }
    SeqFileIn sfi2;
    if (!open(sfi2, toCString( filePath2 ) ) )
    {
        std::cerr << "ERROR: Could not open the second file" << std::endl;
        return 1;
    }

    // read the first sequences from each file object
    StringSet<CharString> metas;
    StringSet<Dna5String> seqHs;
    StringSet<Dna5String> seqVs;
    readRecords(metas, seqHs, sfi1);
    readRecords(metas, seqVs, sfi2);

    Dna5String seqH = concat(seqHs,"",true);
    Dna5String seqV = concat(seqVs,"",true);

    // vector with qgramSizes (for testing)
    std::vector<unsigned> qgramSizes = options.qgramsize;//{31, 15, 10};
    /*
    for (unsigned i = 31; i > 0; --i)
    {
        qgramSizes.push_back(i);
    }*/

    // in this string of seeds we will store the seed chains we will find
    SeedSet<Seed<Simple>> globalSeedSet;

    // define iterator for string of seeds
    typedef Iterator<SeedSet<Seed<Simple> > >::Type TSeedSetIter;

    std::cout << "Length seqH:\t" << length(seqH) << "\tseqV:\t" << length(seqV) << std::endl;

    // we seek for seeds several times. Every next iteration happens in windows where we haven't found any seeds yet.
    // we store those seeds found in repetitive iterations in localSeedChain
    String<Seed<Simple>> localSeedChain;
    Score<int, Simple> scoringScheme( 2,-1,-2 ); // match, mismatch, gap

    std::clock_t start;
    double durationSeeds;
    double durationAlign;
    double durationNW;
    start = std::clock();


    for ( unsigned qGramSize : qgramSizes )
    {
        std::cout << "qgram " << qGramSize << std::endl;
        unsigned i = 0;
        unsigned pos = 0;
        unsigned infixVBegin;
        unsigned infixHBegin;
        unsigned infixVEnd = 0;
        unsigned infixHEnd = 0;
        unsigned nextInfixVBegin = 0;
        unsigned nextInfixHBegin = 0;

        unsigned counter = 0;

        unsigned const seedSetLen = length(globalSeedSet);

        for ( TSeedSetIter globalIter = begin( globalSeedSet, Standard() ); i < seedSetLen + 1; ++i, ++pos )
        {
            clear(localSeedChain);

            infixHBegin = nextInfixHBegin;
            infixVBegin = nextInfixVBegin;

            globalIter = begin( globalSeedSet, Standard() ) + pos;
            if (i == seedSetLen) {
                infixHEnd = length(seqH);
                infixVEnd = length(seqV);
            } else {
                infixHEnd = beginPositionH(*globalIter);
                infixVEnd = beginPositionV(*globalIter);
                nextInfixHBegin = endPositionH(*globalIter);
                nextInfixVBegin = endPositionV(*globalIter);
            }

            //std::cout << i << "/" << seedChainLen << "\tInfix SeqV\t start:\t" << infixVBegin << "\tend:\t" << infixVEnd << std::endl << std::endl;

            //printSeedString(seedChain);

            Dna5String seqVInfix = infix(seqV, infixVBegin, infixVEnd);
            Dna5String seqHInfix = infix(seqH, infixHBegin, infixHEnd);

            if (createSeedChain(localSeedChain,
                                seqHInfix,
                                seqVInfix,
                                qGramSize,
                                Seed<Simple>() ) )
            {
                updateSeedPositions(localSeedChain, infixHBegin, infixVBegin);
                for ( Iterator<String<Seed<Simple> > >::Type localIter = begin( localSeedChain, Standard() ); localIter != end( localSeedChain, Standard() ); ++localIter )
                {
                    if ( !addSeed( globalSeedSet, *localIter, 5 /*max diag dist*/, 10 /*band width*/, scoringScheme, seqH, seqV, Chaos() ) )
                    {
                        addSeed( globalSeedSet, *localIter, Single() ); // just add seed if CHAOS fails
                        ++pos;
                        ++counter;
                    }
                }
            }
        }
    }

    std::cout << "Searching for seeds done" << std::endl;

    durationSeeds = (std::clock() - start )/ (double) CLOCKS_PER_SEC;
    std::cout << "Runtime Find Seeds: " << durationSeeds << std::endl;


    printSeedString(globalSeedSet);

    start = std::clock();

    String<Seed<Simple> > seedChain;
    chainSeedsGlobally( seedChain, globalSeedSet, SparseChaining() );

    std::cout << "global seed chain done" << std::endl;

    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seqH);
    assignSource(row(alignment, 1), seqV);
    AlignConfig<false, false, false, false> alignConfig;  // ordinary global alignment

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    durationAlign = (std::clock() - start )/ (double) CLOCKS_PER_SEC;
    std::cout << "Runtime Banded Chain Alignment: " << durationAlign << std::endl;
    std::cout << "Runtime Total LAGAN: " << durationAlign + durationSeeds << std::endl;

    std::cout << "Score: " << result << std::endl;
    //std::cout << alignment << std::endl;

    //std::cout << toCString(seedChain);

    return 0;
}