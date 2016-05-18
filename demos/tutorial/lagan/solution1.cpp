//![include]
#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>


using namespace seqan;
//![include]

//![createSeedChain]
//![createSeedChainHead]
template <typename TSeedChain, typename TSeqString, typename TQGramSize, typename TSeed>
inline bool createSeedChain(TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, TQGramSize const & qgramSize, TSeed const & /*tag*/)
//![createSeedChainHead]
{
    if ( length(seqH) < qgramSize || length(seqV) < qgramSize ) {
        return false;
    }

    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;
    typedef Infix<Dna5String>::Type TInfix;

    TIndex index(seqH);
    resize(indexShape(index),qgramSize);

    return true;
}
//![createSeedChain]
//![parser]
//![sequences]
struct LaganOptions
{
    CharString path2file1;
    CharString path2file2;
    std::vector<unsigned> qgramsize;
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
//![parser]
//![main]
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

    //vector for qgram sizes
    std::vector<unsigned> qgramSizes = options.qgramsize;
    //![sequences]
    return 0;
}

//![main]