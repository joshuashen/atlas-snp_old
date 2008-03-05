## 

fa = ARGV[0]
ss = ARGV[1].to_i
pp = ARGV[2]


flags = {}

files = []

num = 1

fo = pp+'_'+num.to_s+ '.fa'
fout = File.new(fo, "w")
files << fo

seq = 0
File.new(fa,'r').each do |line|
  if line=~ /^>(\S+)/
    name = $1
    if seq >= ss
      flags[name] = 1
      seq = 0
      fout.close
      num += 1
      fo = pp+'_'+num.to_s+ '.fa'
      fout = File.new(fo, "w")
      files << fo
    end
    
    fout.puts line
  else
    fout.puts line
    seq += line.chomp.size
  end
end

fout.close

puts files.join("\n")
