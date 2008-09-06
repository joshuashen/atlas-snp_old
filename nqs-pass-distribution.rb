#!/usr/bin/env ruby
require 'getoptlong'
require 'bigdecimal'

opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

homeDir=Dir.pwd
outputDir=homeDir+'/output/'

file_name=File.basename(optHash["--input"])
file_path=File.dirname(optHash["--input"])

file_output=outputDir+file_name+'.nqs'
fo_output=File.new(file_output,'w')
#$stderr.puts file_output

file_dis=outputDir+file_name+'.nqs.dis'
fo_dis=File.new(file_dis,'w')
#$stderr.puts file_dis

$countofReads=0
$count=0
$nqs={}
name=nil
$qual=[]
$nqs= Hash.new {|h,k| h[k]=Hash.new}
File.new(optHash["--input"],'r').each do |line|
  if line=~ /^>(\S+)/
    if $countofReads>0
      $nqs[name][:totalCount]=$count
      $nqs[name][:elem]=$qual
      $nqs[name][:nqs][0]="1"
      $nqs[name][:nqs][1]="1"
      $nqs[name][:nqs][2]="1"
      $nqs[name][:nqs][3]="1"
      $nqs[name][:nqs][4]="1"
      $nqs[name][:passed]+=5
      $nqs[name][:nqs][$count-1]="NA"
      $nqs[name][:nqs][$count-2]="NA"
      $nqs[name][:nqs][$count-3]="NA"
      $nqs[name][:nqs][$count-4]="NA"
      $nqs[name][:nqs][$count-5]="NA"
      if $nqs[name][:totalCount]>10
        i=5
        ($nqs[name][:totalCount]-10).times do
          if $qual[i-1].to_i>=15 and $qual[i-2].to_i>=15 and $qual[i-3].to_i>=15 and $qual[i-4].to_i>=15 and $qual[i-5].to_i>=15 and $qual[i+1].to_i>=15 and $qual[i+2].to_i>=15 and $qual[i+3].to_i>=15 and $qual[i+4].to_i>=15 and $qual[i+5].to_i>=15 and $qual[i].to_i>=20
            $nqs[name][:passed]+=1
          else
            $nqs[name][:failed]+=1
            $nqs[name][:nqs][i]=0
          end
          i+=1
        end
        $nqs[name][:ratio]=($nqs[name][:passed].to_f/($nqs[name][:totalCount]-5).to_f).to_f
      else
        #$nqs[name][:ratio]=$qual
        j=0
        sum=0
        ($nqs[name][:totalCount]).times do
          sum+=$qual[j]
        end
            $nqs[name][:nqs][i]=1
        $nqs[name][:ratio]=sum.to_f/$nqs[name][:totalCount].to_f
      end
      #print "#{name}\t#{$nqs[name][:passed]}\t#{$nqs[name][:totalCount]}\t#{$nqs[name][:ratio]}\n"
    end
    
    $count=0
    name=$1
    $nqs[name][:totalCount]=0
    $nqs[name][:passed]=0
    $nqs[name][:failed]=0
    $nqs[name][:ratio]=0
    $nqs[name][:nqs]={}
    $countofReads+=1
    $qual=[]
  else
   cols=line.split(/\s+/)
   $count+=cols.length
   cols.each do |elem|
    $qual << elem
   end  
  end
end
$nqs[name][:totalCount]=$count
$nqs[name][:elem]=$qual
$nqs[name][:nqs][0]="1"
$nqs[name][:nqs][1]="1"
$nqs[name][:nqs][2]="1"
$nqs[name][:nqs][3]="1"
$nqs[name][:nqs][4]="1"
$nqs[name][:passed]+=5
$nqs[name][:nqs][$count-1]="NA"
$nqs[name][:nqs][$count-2]="NA"
$nqs[name][:nqs][$count-3]="NA"
$nqs[name][:nqs][$count-4]="NA"
$nqs[name][:nqs][$count-5]="NA"

if $nqs[name][:totalCount]>10
  i=5
  ($nqs[name][:totalCount]-10).times do
    if $qual[i-1].to_i>=15 and $qual[i-2].to_i>=15 and $qual[i-3].to_i>=15 and $qual[i-4].to_i>=15 and $qual[i-5].to_i>=15 and $qual[i+1].to_i>=15 and $qual[i+2].to_i>=15 and $qual[i+3].to_i>=15 and $qual[i+4].to_i>=15 and $qual[i+5].to_i>=15 and $qual[i].to_i>=20
      $nqs[name][:passed]+=1
    else
      $nqs[name][:failed]+=1
      $nqs[name][:nqs][i]=0
    end
    i+=1
  end
  $nqs[name][:ratio]=($nqs[name][:passed].to_f/($nqs[name][:totalCount]-5).to_f).to_f
else   
  j=0
  sum=0
  ($nqs[name][:totalCount]).times do
    sum+=$qual[j]
  end
  $nqs[name][:nqs][i]=1
  $nqs[name][:ratio]=sum.to_f/$nqs[name][:totalCount].to_f
end
    
$nqs.keys.sort.each do |name|
  fo_output.print "#{name}\t#{$nqs[name][:passed]}\t#{$nqs[name][:totalCount]}\t#{$nqs[name][:ratio]}\n"
end

fo_output.close
$stderr.puts 'nqs for'+file_name+' has been done!'
#####Distribution#############################
fo_read_output=File.new(file_output,'r')
$count=[]

0.upto(9) do |i|
  $count[i]=[]
  $count[i]<<0
  $count[i]<<0
  $count[i]<<0
  $count[i]<<0
end

#$stderr.puts $count.size
zero=0
hash_passRatio={}

fo_read_output.each do |line|
  cols=line.split(/\s+/)
  
  query=cols[0]
  passRatio=cols[3].to_f
  
  hash_passRatio[query]=passRatio
  if passRatio==0
    zero+=1
  end
  
  0.upto(9) do |i|
    0.upto(3) do |j|
      j1=j+1
      v1=i.to_f/10+j1.to_f/40
   
      p=BigDecimal(passRatio.to_s)
      v=BigDecimal(v1.to_s)
      if p==v
        $count[i][j]+=1
      end
    end
  end
end

###Print result###########
#$stderr.puts count[1]1
$counter=[]
0.upto(9) do |i|
  $counter[i]=[]
  $counter[i]<<0
  $counter[i]<<0
  $counter[i]<<0
  $counter[i]<<0
end

0.upto(9) do |i|
  0.upto(3) do |j|
    half=$count[i][j]/2
    #$stderr.puts half
    if $count[i][j]%2==0
      if j<3
        $counter[i][j]+=half
        $counter[i][j+1]+=half
      elsif j==3 and i<9
        $counter[i][j]+=half
        $counter[i+1][j-3]+=half
      elsif j==3 and i==9
        $counter[9][3]+=$count[i][j].to_i
      end
    else
      if j<3
        $counter[i][j]+=half
        $counter[i][j+1]+=half+1
      elsif j==3 and i<9
        $counter[i][j]+=half
        $counter[i+1][j-3]+=half+1
      elsif j==3 and i==9
        $counter[9][3]+=$count[i][j].to_i
      end
    end
  end
end

hash_passRatio.keys.sort.each do |q|
  ratio=hash_passRatio[q]
  0.upto(9) do |i|
    0.upto(3) do |j|
      v1=i.to_f/10+j.to_f/40
      v2=i.to_f/10+j.to_f/40+0.025
      r=BigDecimal(ratio.to_s)
      v11=BigDecimal(v1.to_s)
      v22=BigDecimal(v2.to_s)
      #$stderr.puts v1.to_s+" "+v2.to_s
      if r>v11 and r<v22
        $counter[i][j]+=1
      end
    end
  end

end

c=0
$counter[0][0]+=zero
while $counter.size>0
  a=$counter.shift
  a1=a.shift
  a2=a.shift
  a3=a.shift
  a4=a.shift
  
  fo_dis.print "#{c.to_f/10}-#{c.to_f/10+0.025}\t#{a1}\n"
  fo_dis.print "#{c.to_f/10+0.025}-#{c.to_f/10+0.05}\t#{a2}\n"
  fo_dis.print "#{c.to_f/10+0.05}-#{c.to_f/10+0.075}\t#{a3}\n"
  fo_dis.print "#{c.to_f/10+0.075}-#{c.to_f/10+0.1}\t#{a4}\n"
  c+=1
end

fo_dis.close

$stderr.puts "NQS distribution has been done!"
####write to file
#$stderr.puts $countofReads


