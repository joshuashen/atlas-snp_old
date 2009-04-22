#!/usr/bin/env ruby
require 'getoptlong'
require 'bigdecimal'

#Program overview:
#       'Nqs.pass.rb' is used to determine whether a variant base pass Neighborhood Quality Standard (NQS).
#   A base is considered as passing the NQS window threshold if both its own “phred-like” quality score is
#   at least 20, and the quality scores of each of the five flanking bases on either side are at least 15.
#   This program takes three parameters. The first is the fasta file containing quality scores of all 454
#   bases introduced by an '--qual'option.  The second parameter is the crossmatch file introduced by
#   an '--xm' option. The option '--outputPrefix' is defined to let user to set the prefix name of the output.
#   When running this program, an output file
#   named 'outputPrefix.NQS' containing the 'NQS pass' information for each variant base will be created.
#
#Command Line Arguments:
#  Usage: nqsPass.rb -i [fasta file] -d [file.dist] -o [outputfile root directory]
#  -q is a required argument, which must be followed by fasta file of quality score.
#  -x is a required argument, associated with the crossmatch file.
#  -o is a required argument which must be appended by the user-defined output root directory without the last slash.
#
#Example:
#   Create NQS file for quality fasta file 'fasta.part_0_1.fa.qual' by taking crossmatch file 'chr1.xm.part_0'
#  as the input. The output file root directory will be 'chr1'.
#       ruby nqsPass_op_V3.1.rb -i fasta.part_0_1.fa.qual -d reads.part_0.dist -o chr1
#
#Output format:
#  'read_name[tab]read_length[tab]distance_to_3'[tab]quality[tab]NQS_pass
class  NqsPass

#Initialize options 
def NqsPass.processArguments()

opts = GetoptLong.new(
    ["--qual", "-q", GetoptLong::REQUIRED_ARGUMENT],
    ["--xm", "-x", GetoptLong::REQUIRED_ARGUMENT],
    ["--help","-h", GetoptLong::NO_ARGUMENT],
    ["--outputPrefix","-o", GetoptLong::REQUIRED_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

NqsPass.usage() if (optHash.key?("--help"));

NqsPass.usage() if (optHash.empty?);
return optHash
end

#The usage information of the program
def NqsPass.usage(msg='')
    unless (msg.empty?)
        puts "\n#{msg}\n"
    end
    puts "
Program overview:
       'Nqs.pass.rb' is used to determine whether a variant base pass Neighborhood Quality Standard (NQS).
   A base is considered as passing the NQS window threshold if both its own “phred-like” quality score is
   at least 20, and the quality scores of each of the five flanking bases on either side are at least 15.
   This program takes three parameters. The first is the fasta file containing quality scores of all 454 bases introduced
   by an '--qual'option.  The second parameter is the crossmatch file introduced by an '--xm' option.
   The option '--outputFileRoot' is defined to let user denote their working directory without the last '/'.
   Thus, the program will automatically create an output directory '_.nqs-output' under the working directory as shown in the argument.
   After running this program, an output file with suffix '.NQS' containing the 'NQS pass' information for each
   variant base will be created under the output directory.

Command Line Arguments:
  Usage: nqsPass.rb -i [fasta file] -d [file.dist] -o [outputfile root directory]
  -i is a required argument, which must be followed by fasta file of quality score.
  -x is a required argument, associated with the crossmatch file.
  -o is a required argument which must be appended by the user-defined output root directory without the last slash.
  
Example:
  Create NQS file for quality fasta file 'fasta.part_0_1.fa.qual' by taking crossmatch file 'chr1.xm.part_0'
  as the input. The output file root directory will be 'chr1'.
    ruby nqsPass_op_V3.1.rb -i fasta.part_0_1.fa.qual -x chr1.xm.part_0 -o chr1

Output format:
  'read_name[tab]read_length[tab]distance_to_3'[tab]quality[tab]NQS_pass"
    exit(2);
end

#Initialize the object of the class "NqsPass"
def initialize(optHash)
    @optHash=optHash
    setParameters()
end

#Initialize the parameters
def setParameters()
    if (!@optHash.key?('--qual') or !@optHash.key?("--xm") or !@optHash.key?("--outputPrefix") )then
        NqsPass.usage("Option missing!")
        exit(2);
    end
    
    @xmFile=@optHash["--xm"]
    @fastaFile=@optHash["--qual"]
    #@outputRoot=@optHash["--outputFileRoot"]
    @outputPrefix=@optHash["--outputPrefix"]
    
    @file_name=File.basename(@fastaFile)
    @file_path=File.dirname(@fastaFile)

    #Initialize the output directory
    #@outputDir=@outputRoot+'/'+@file_name+'.nqs-output/'
    #system("mkdir #{@outputDir}") #Create the directory

    #Initialize the files containing NQS information
    @file_NQS=@outputPrefix+'.NQS'
end

$pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/

#Parse out the variation position from crossmatch file
def extractVariationPosition()
    @hash_fasta_dist={}
    query,qsize,score= '',0, 0
    sub= 0
    xmReader=File.open(@xmFile,"r")
    xmReader.each do |line|
    if line.match($pattern)
        score,sub,del,ins,query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14
        qsize = qend + qright
   
        elsif line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\((\d+)\)\s+(\d+)\s+(\S+)/
          type, qplace, tplace = $1, $2.to_i, $5.to_i
    
          if type =~ /^S/ # substitution
            dist = qsize - qplace  #calcuate distance to 3'
            if @hash_fasta_dist[query]==nil
                @hash_fasta_dist[query]=[]
            end
            @hash_fasta_dist[query]<<dist 
          end
        end
    end
    $stderr.puts @hash_fasta_dist.size
    xmReader.close
end

#Determine 'NQS pass' status.
def calculateNQSPass()
    @hash_nqs={}
    p_status=''
    name='1'
    qual=[]
    count=0
   fastaReader=File.open(@fastaFile,"r")
   puts @fastaFile
   
   fastaReader.each do |line|
     if line=~ /^>(\S+)/
        if @hash_fasta_dist[name]!=nil
                @hash_fasta_dist[name].uniq.sort.each do |i|
                    dist=i.to_i
                    #puts dist
                    pos=count-1-dist #Calculate the base position
                    if count>10 #Read length >= 11
                        if dist <5 #Last four bases
                            p_status="0" #Set NQS pass to be "0", failed
                        elsif pos<5 #The first four bases
                            p_status="1" #Set NQS pass to be "1", passed
                        else #Considerring the bases having five flanking bases
                            #Judge whether NQS passes
                            #puts qual[pos].to_s+"###"
                            if qual[pos].to_i>=20 and qual[pos-1].to_i>=15 and qual[pos-2].to_i>=15 and qual[pos-3].to_i>=15 and qual[pos-4].to_i>=15 and qual[pos-5].to_i>=15 and qual[pos+1].to_i>=15 and qual[pos+2].to_i>=15 and qual[pos+3].to_i>=15 and qual[pos+4].to_i>=15 and qual[pos+5].to_i>=15
                                
                                p_status="1" #Set NQS pass to be "1", passed
                            else 
                                p_status="0" #Set NQS pass to be "0", failed
                            end
                        end
                    else #Read length <= 10
                        p_status="1" #Set NQS pass to be "1", passed
                    end
                    
                    #Print the distance to 3' (pos), "NQS pass", Quality and Read length into a string
                    str=dist.to_s+'.'+p_status+'.'+qual[pos].to_s+'.'+count.to_s
                
                    @hash_nqs[name]<<str
                end
           
        end
        ##Initialize the hash table for storing NQS information
        name=$1 #Read name
        next if @hash_fasta_dist[name]==nil #Ignore those reads not existing in .dist file
        if @hash_fasta_dist[name]!=nil 
            if @hash_nqs[name]==nil
                @hash_nqs[name]=[]
            end
        end
        count=0 #Counter of read length
        qual=[] #array of quality scores
    else
        next if @hash_fasta_dist[name]==nil
        if @hash_fasta_dist[name]!=nil
            #start to write to quality array
            cols=line.split(/\s+/)
            count+=cols.length #Adding to counter
            #cols.each do |elem|
            #  $qual << elem
            #  
            #end
            qual.concat(cols)
        end
    end  
   end
   if @hash_fasta_dist[name]!=nil
             #   $stderr.puts hash_fasta_dist[name].size
                @hash_fasta_dist[name].uniq.sort.each do |i|
                    dist=i.to_i
                    pos=count-1-dist
                    if count>10
                        if dist <5
                            p_status="0"
                        elsif pos<5
                            p_status="1"
                        else
                            pos_center=qual[pos].to_i
                            if qual[pos].to_i>=20 and qual[pos-1].to_i>=15 and qual[pos-2].to_i>=15 and qual[pos-3].to_i>=15 and qual[pos-4].to_i>=15 and qual[pos-5].to_i>=15 and qual[pos+1].to_i>=15 and qual[pos+2].to_i>=15 and qual[pos+3].to_i>=15 and qual[pos+4].to_i>=15 and qual[pos+5].to_i>=15
                            
                                p_status="1"
                            else
                                p_status="0"
                            end
                        end
                    else
                        p_status="1"
                    end
                    str=dist.to_s+'.'+p_status+'.'+qual[pos].to_s+'.'+count.to_s
                    
                    @hash_nqs[name]<<str
                end
            #end
        end
   $stderr.puts @hash_nqs.size
   fastaReader.close
end

#Print result
def printResult()
    nqsWriter=File.new(@file_NQS,"w")
    nqsWriter.print "readName\treadLength\tdistanceToRightEnd\tQuality\tnqsPASS\n"
    @hash_nqs.keys.sort.each do |n|
        @hash_nqs[n].sort.each do |i|
            str=i
            if str=~ /(\S+).(\S+)\.(\S+)\.(\S+)/
              d=$1
              p=$2
              qual=$3
              length=$4
            end
            nqsWriter.print "#{n}\t#{length}\t#{d}\t#{qual}\t#{p}\n"
        end
    end
    nqsWriter.close
end
end

#Main methods
optHash=NqsPass.processArguments()
nqsPassDeterminter=NqsPass.new(optHash)
nqsPassDeterminter.extractVariationPosition()
nqsPassDeterminter.calculateNQSPass()
nqsPassDeterminter.printResult()
exit(0);