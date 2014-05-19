Enter file contents here
#!/usr/bin/perl             
$aa=shift @ARGV;             #该程序用于从fasta文件中根据指定的位置信息，抓取序列；第一个参数是位置信息，每一行格式如下：contig_35424:2:207（用冒号分开，序列ID:位置1：位置2）
$bb=shift @ARGV;             #第二个参数是fasta文件
open AA,"<$aa";               
$aa=~s/(\w*)\.txt/$1.fasta/;            #$aa作为输出文件的文件名
open CC,">>$aa";
while(<AA>){                       
  chomp;
  @array=split /:/,$_;                  #根据冒号风格字符串到数组中
  $flag=1;
  open BB,"<$bb";
  while(chomp($string=<BB>)){           #根据数组中的信息，从打开的fasta搜索序列，并保存
   if($flag&&$string=~/$array[0]/){     #变量flag初始为1，便于寻找序列ID(eg：contig_35424)
         $flag=0;                       #当寻找到序列ID后，flag设为0,之后开始寻找特定字符串
	     next;                          #类似c与pathon中的continue； 
     }
   elsif(!$flag&&$string!~/>.*/){        
         $temp_string=length($string);   #判断字符串长度
         if($array[1]>$temp_string){     #如果位置1坐标大于字符串长度，那么将位置减去当前字符串长度，并读入下一行
             $array[1]-=$temp_string;  
			 $array[2]-=$temp_string;
		 }
		 elsif($array[2]>=$temp_string){       #位置2坐标大于当前字符串长度，则开始读取子字符串
		     if($array[1]>0){
		     $mRNA.=substr($string,$array[1]-1);      #在字符串中寻找字串
			 $array[1]-=$temp_string;                 #寻找到匹配的字符串后，将位置的坐标相应将去当前字符串长度
			 $array[2]-=$temp_string;
		    }
			else{ 
			 $mRNA.=substr($string,0);
			 $array[2]-=$temp_string;
			 }
		 }
		 elsif($array[2]>0){ 
        	 if($array[1]>0){ 
			     $mRNA.=substr($string,$array[1]-1,$array[2]-$array[1]+1);
		         $array[2]-=$temp_string;
		     }
			 else{  
			     $mRNA.=substr($string,0,$array[2]);
		         $array[2]-=$temp_string;
			 }
		 }
		 if($array[2]>0){next;}
	}
	else{ next;}
	if($array[2]>0){ print CC ">$array[0]\nthis string no exists\n";last;}
	print CC ">$array[0]\n$mRNA\n";
	$mRNA='';                                    #将mRNA重新设为空；
	last;
  }
  close BB;                                      #关闭文件句柄
}
close AA;
close CC;
