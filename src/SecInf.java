class SecInf
{  
	public double u;                    // u-value
   	public Solid solid;                 // primitive solid hit
   	public SecInf next;                 // link
   	
   	public SecInf(double u, Solid s, SecInf n)
   	{
   		this.u = u;
   		solid = s;
   		next = n;
   	}
}
