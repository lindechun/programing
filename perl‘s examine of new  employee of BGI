#!/usr/bin/perl
#use Time::HiRes qw/time/;       #计算运行时间的模块
#$t_start=time();
unless(-f jieguo){mkdir "./jieguo",0755 or warn "can't make fred directory:$!" };  #创建目录用于存放结果文本
print "the first question that calculate those length begin\n------------Press Enter to start------------\n";
$a=<>;

#第一题：计算长度并输出到文本
open IN,"<plum_0630.scafSeq.FG" or die("Can't find file",$!);
while (<IN>){
    if($.%2==0){                                                          #对文件中的偶数行统计序列信息
	     $TotalLen+=tr/ATCGN/ATCGN/;                                      #返回替换次数
	     $EffectiveLen+=tr/ATCG/ATCG/;
	     $NLen+=tr/N/N/;
	     $GCLen+=tr/GC/GC/;
    }
}
$GCRate=100*$GCLen/$EffectiveLen;
$GCRate=sprintf "%0.2f",$GCRate;
open OUT, ">./jieguo/plum_0630.scafSeq.txt" or die("Can't write file",$!);             #输出到文本
print OUT "Name\tTotal Length\tEffective Length\tN Length\tGC Length\tGC Rate(%)\t\n";
print OUT "scaffold.fa\t$TotalLen\t$EffectiveLen\t$NLen\t$GCLen\t$GCRate\n";
close IN,OUT;
print "the second question that calculate N50,N90 begin\n------------Press Enter to start------------\n";
$a=<>;



##第二题：计算N50、N90 scaffold,输出到文本
open IN,"<plum_0630.scafSeq.FG" or die("Can't find file",$!);
open OUT,">./jieguo/N50_N90.txt" or die("Can't find file",$!);
my @scaffold;                                #保存到数组@scaffold和@sequence
my @sequence;
while(chomp(my $line=<IN>)){
	push @scaffold,$line;
    chomp($line=<IN>);
    push @sequence,$line;
}
my @LenArr;                                  #计算Sequence长度,保存到数组@LenArr
my %hashscaff;
foreach my $e(@sequence){
	 push @LenArr,length $e;        ###多行注释处是原始版，此处修改过
     $hashscaff{length $e}=$e;
=plob
     $Len=$e=~tr/ATCGN/ATCGN/;
	 push @LenArr,$Len;
     $hashscaff{$Len}=$e;
=cut
}
my @sortedArr=sort {$b<=>$a} @LenArr;           #按长度从大到小排序，保存到数组@sortedArr
$N50Len=$TotalLen*0.5;                          #计算总长的50%和90%,保存到$N50Len,$N90Len
$N90Len=$TotalLen*0.9;
for($i=0,my $Len50=0;$Len50<$N50Len;$i++){      #累加求取N50scaffold,输出到文本
	$Len50+=$sortedArr[$i];
}
print OUT "N50 scaffold is $sortedArr[$i]bp\n";
print OUT "$hashscaff{$sortedArr[$i]}\n";                          
for($i=0,my $Len90=0;$Len90<$N90Len;$i++){       #累加求取N90scaffold，输出到文本
	$Len90+=$sortedArr[$i];
}
print OUT "N90 scaffold is $sortedArr[$i]bp\n";
print OUT "$hashscaff{$sortedArr[$i]}\n";
close IN,OUT;
print "the third question that calculate GC  percent  begin\n------------Press Enter to start------------\n";
$a=<>;



##第三题：以250bp为窗口，无重叠滑动，GC含量
open OUT,">./jieguo/GCcontent.txt" or die $!;
for($i=0;$sequence[$i];$i++){
    print OUT $scaffold[$i]."\n";
	for($n=0;$n<$LenArr[$i];$n+=250){
        $window=substr($sequence[$i],$n,250);
        $GCLen=$window=~tr/GC/GC/;
        $windowLen=$window=~tr/ATCG/ATCG/;
        if($windowLen ne 0){
            $GCcontent=$GCLen/$windowLen;  
            printf OUT "%.4f\t",$GCcontent;
        }
	}
    print OUT "\n";
}
close OUT;
print "the fourth question that search CDS sequence  begin\n------------Press Enter to start------------\n";
$a=<>;



##第四题：输出CDS到fa文件
open IN,"<Prunus_mume_scaffold.gff" or die$!;
open OUT,">./jieguo/Prunus_mume_scaffold.fa" or die $!;
my $strand;
my $seq;
my $flag=0;
while(my $line=<IN>){                   #此处的while与后面多行注释掉的while相比，运行时间更短，因为$TempArr[0]变量匹配到数组元素后,
	@TempArr=split/\t/,$line;           #使用$flag记住状态，这样遇到下一个CDS时，直接在对应的数组元素里取出字符串
	if($TempArr[2] eq "CDS"){
       if(!$flag){
		for($i=0;$scaffold[$i];$i++){
			if($scaffold[$i]=~/$TempArr[0] /){     #此处的模式匹配最后的空格不能删除，因为scaffold16与scaffold162是不同的
                $flag=1;$k=$i;
			    $CDS=substr($sequence[$i],$TempArr[3]-1,$TempArr[4]-$TempArr[3]+1);         #取出CDS
			    $seq=$seq.$CDS;
                last;
			}
		}
       }else{   
           $CDS=substr($sequence[$k],$TempArr[3]-1,$TempArr[4]-$TempArr[3]+1);
           $seq=$seq.$CDS;
       }
	}else{
		if($seq ne ""){
             if($strand eq "-"){$seq=reverse $seq;$seq=~tr/ATCG/TAGC/;}
		print OUT "$seq\n";
		$seq="";
        }
	    print OUT ">$TempArr[0]\t$TempArr[8]";
        $strand=$TempArr[6];
        $flag=0;
    }
}
=plod                               ###########多行注释
while(my $line=<IN>){
	@TempArr=split/\t/,$line;
	if($TempArr[2] eq "CDS"){
		for($i=0;$scaffold[$i];$i++){
			if($scaffold[$i]=~/$TempArr[0] /){
			    $CDS=substr($sequence[$i],$TempArr[3]-1,$TempArr[4]-$TempArr[3]+1);         #取出CDS
			    $seq=$seq.$CDS;
                last;
			}
		}
	}else{
		if($seq ne ""){
             if($strand eq "-"){$seq=reverse $seq;$seq=~tr/ATCG/TAGC/;}
		print OUT "$seq\n";
		$seq="";
        }
	    print OUT ">$TempArr[0]\t$TempArr[8]";
        $strand=$TempArr[6];
    }
}
=cut
if($strand eq "-"){$seq=reverse $seq;$seq=~tr/ATCG/TAGC/;}
print OUT "$seq\n";
close IN,OUT;
print "the fifth  question that translate amino acid begin\n------------Press Enter to start------------\n";
$a=<>;



##第五题：翻译CDS，输出氨基酸序列到文件
my %aminoHash = (    #构建密码子的哈希表
    'ACG' => 'T',    # Threonine   
    'ACT' => 'T',    # Threonine   
    'AAC' => 'N',    # Asparagine   
    'AAT' => 'N',    # Asparagine   
    'AAA' => 'K',    # Lysine   
    'AAG' => 'K',    # Lysine   
    'AGC' => 'S',    # Serine   
    'AGT' => 'S',    # Serine   
    'AGA' => 'R',    # Arginine   
    'AGG' => 'R',    # Arginine   
    'GTA' => 'V',    # Valine   
    'GTC' => 'V',    # Valine   
    'GTG' => 'V',    # Valine   
    'GTT' => 'V',    # Valine   
    'GCA' => 'A',    # Alanine   
    'GCC' => 'A',    # Alanine   
    'GCG' => 'A',    # Alanine   
    'GCT' => 'A',    # Alanine       
    'GAC' => 'D',    # Aspartic Acid   
    'GAT' => 'D',    # Aspartic Acid   
    'GAA' => 'E',    # Glutamic Acid   
    'GAG' => 'E',    # Glutamic Acid   
    'GGA' => 'G',    # Glycine   
    'GGC' => 'G',    # Glycine   
    'GGG' => 'G',    # Glycine   
    'GGT' => 'G',    # Glycine      
);
open IN,"<./jieguo/Prunus_mume_scaffold.fa" or die("Can't find file",$!);
open OUT,">./jieguo/Amino_acid.fa" or die $!;
while(chomp(my $line=<IN>)){
	if($line=~/scaffold/){
         print OUT "$line\n"
    }else{
		$len=length $line;
		for(my $s=0;$s<$len;$s+=3){
            print OUT $aminoHash{substr($line,$s,3)};   #将下面被注释了的两行合并成一行，
            #$code=substr($line,$s,3);       
			#print OUT $aminoHash{$code};
		}
		print OUT "\n";
	}
}
close IN,OUT;
print "Amino_acid.fa established\n\nall question have solve\n------------Press any key to end------------\n";
$a=<>;
#my $t_end   = time();
#open TIME ,">time.txt";
#printf TIME "used: %.6f sec\n", $t_end - $t_start;
