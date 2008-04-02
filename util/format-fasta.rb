## limit the length of each line to 80

input = ARGV[0]

output = input + '.mod.fa'

$fin = File.new(input, 'r')

$fout = File.new(output, 'w')


def format(seq)
  return if seq == ''
  s = seq.size
  0.upto(s/80) do |i|
    $fout.syswrite(seq[i*80..(i+1)*80-1] + "\n")
#    $fout.syswrite(seq.slice!(0,80) + "\n")
  end
end

str = ''

$fin.each do |line|
  if line=~ /^>(\S+)/
    format(str)
    $fout.syswrite(line)
    str = ''
  else
    str << line.chomp
  end
end

format(str)

$fin.close
$fout.close

