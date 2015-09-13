difference ()
{
	import("N27E086_join_trim.stl",convexity=10);
	linear_extrude(height=97.5,convexity=10)
	{
		union()
		{
			offset(r=0.5)
			{
				import("Everest_1953_Route.dxf",convexity=10);
			}
			offset(r=-0.5)
			{
				import("Everest_1953_Route.dxf",convexity=10);
			}
		}
	}
}
