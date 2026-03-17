#/usr/bin/perl
use strict;

# 定义哈希存储 SNPs
my %snps_hash;

# 获取参数
my $snps = $ARGV[0] or die "Usage: $0 <snp_list> <tissue> <annot_file>\n";
my $tissue = $ARGV[1] or die "Usage: $0 <snp_list> <tissue> <annot_file>\n";
my $annot_file = $ARGV[2] or die "Usage: $0 <snp_list> <tissue> <annot_file>\n";

# 检查文件是否存在
die "SNP list file not found: $snps" unless -e $snps;
die "Annotation file not found: $annot_file" unless -e $annot_file;

# 读取 SNP 列表
open my $snps_file, '<', $snps or die "Failed to open SNP list: $!\n";
while (<$snps_file>) {
  chomp;
  $snps_hash{$_} = 1;
}
close $snps_file;

# 读取注释文件
open my $annot_fh, '<', $annot_file or die "Failed to open annotation file: $!\n";

# 打印表头
my $header = <$annot_fh>;
chomp $header;
print "CHR\tBP\tSNP\tCM\t${tissue}\n";

# 检查表头格式
my @header_cols = split '\t', $header;
die "Annotation file header format error: $header" unless @header_cols >= 4;

# 处理每一行
while (<$annot_fh>) {
  chomp;
  my @cols = split '\t';
  die "Annotation file row format error: $_" unless @cols >= 4;
  my ($chr, $bp, $snp, $cm) = @cols[0..3];
  print "$chr\t$bp\t$snp\t$cm\t", exists($snps_hash{$snp}) ? 1 : 0, "\n";
}
close $annot_fh;