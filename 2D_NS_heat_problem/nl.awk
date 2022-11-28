#
#	Have to add a newline for a new row of coordinates
#
BEGIN { OFS="	"; YO=-1.23456789; X=YO; Y=YO; Z=YO }
      { 
        if ($1!="")
	{
	  if ($1!=YO) { print " "; YO=$1 }
          if ($1==X && $2==Y)
	   {
#	     print $1,$2,($3+Z)/2
	   }
	   else
	   {
	     print $1,$2,$3
	   }
	  X=$1; Y=$2; Z=$3;
	}
      }
END {}
