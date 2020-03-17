package roadgraph;

import java.util.Comparator;

public class AStarDistComparator implements Comparator<MapNode> {
	public int compare(MapNode n1, MapNode n2) {
		if (n1.getAStarDistance() < n2.getAStarDistance()) {
			return -1;
		}
		else if (n1.getAStarDistance() > n2.getAStarDistance()) {
			return 1;
		}
		return 0;
	}
}
