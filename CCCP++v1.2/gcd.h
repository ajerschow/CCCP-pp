// Euclid's algorithm

using namespace std;

int gcd( int num1, int num2 )
{
  if (num1==0)
    return num2;
  else {
    int remainder = num2 % num1;
    
    if ( remainder != 0 )
      return gcd( remainder,num1 );
    
    return num1;
  }
}
                                            
int lcm( int num1, int num2)
{
  if (num1==0 || num2==0)
    return 0;
  else
    return num1*num2/gcd(num1,num2);
}



////////////////////////////////////
// CAREFUL ABOUT GCD ROUTINES     //
////////////////////////////////////
int gcd(const vector<int>& nums)
{
  // CAREFUL: IF SIZE=1 RETURNS the number

  int i,result=nums[0];
  if(nums.size()==1)
    return nums[0];
  else
    for (i=1;i<nums.size();i++)
      result=gcd(result,nums[i]);
  
  return result;
}

int lcm(const vector<int>& nums)
{
  int i,result=nums[0];
  if(nums.size()==1)
    return nums[0];
  else
    for (i=1;i<nums.size();i++)
      result=lcm(result,nums[i]);
  return result;
}

//gcd(m1,m2)lcm(m1,m2) = m1m2
//gcd(m1, (gcd(m2, m3)) = gcd(gcd(m1, m2), m3) 
//lcm(m1, (lcm(m2, m3)) = lcm(lcm(m1, m2), m3) 
