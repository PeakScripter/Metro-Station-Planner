package jminor;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.util.*;
import java.util.List;
import java.awt.geom.Point2D;
class MetroStationPlannerGUI {
    private jminor.MetroStationPlanner planner;
    private MetroClusterVisualizer visualizer;
    private Map<Integer, Point2D> clusterCenters;
    private Set<Integer> finalStations;

    public MetroStationPlannerGUI() {
        planner = new jminor.MetroStationPlanner();
        clusterCenters = new HashMap<>();
        finalStations = new HashSet<>();

        // Set up main frame
        JFrame frame = new JFrame("Metro Station Planner");
        frame.setSize(1000, 700);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(null);

        // File selection button
        JButton fileButton = new JButton("Load Dataset");
        fileButton.setBounds(20, 20, 150, 30);

        // Preprocess button
        JButton preprocessButton = new JButton("Preprocess Data");
        preprocessButton.setBounds(20, 60, 150, 30);
        preprocessButton.setEnabled(false);

        // Cluster button
        JButton clusterButton = new JButton("Apply DBSCAN");
        clusterButton.setBounds(20, 100, 150, 30);
        clusterButton.setEnabled(false);

        // Betweenness Centrality button
        JButton centralityButton = new JButton("Calculate Stations");
        centralityButton.setBounds(20, 140, 150, 30);
        centralityButton.setEnabled(false);

        // Panel for graph visualization
        JPanel graphPanel = new JPanel();
        graphPanel.setBounds(200, 20, 760, 620);
        graphPanel.setLayout(new BorderLayout());
        graphPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK));

        // Action listener for file selection
        fileButton.addActionListener(e -> {
            JFileChooser fileChooser = new JFileChooser();
            if (fileChooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
                try {
                    planner.loadData(fileChooser.getSelectedFile().getAbsolutePath());
                    JOptionPane.showMessageDialog(frame, "Data Loaded Successfully");
                    preprocessButton.setEnabled(true);
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(frame, "Error loading file: " + ex.getMessage());
                }
            }
        });

        // Action listener for preprocessing
        preprocessButton.addActionListener(e -> {
            try {
                planner.preprocessData(planner.getData());
                JOptionPane.showMessageDialog(frame, "Data Preprocessed Successfully");
                clusterButton.setEnabled(true);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(frame, "Error in preprocessing: " + ex.getMessage());
            }
        });

        // Action listener for clustering
        clusterButton.addActionListener(e -> {
            try {
                double[][] distanceMatrix = planner.customGeographicalDistanceMatrix(planner.getData());
                int[] clusters = planner.applyDBSCAN(distanceMatrix, 125, 2);
                clusterCenters = calculateClusterCenters(planner.getData(), clusters);

                // Calculate MST edges
                Set<Line2D> mstEdges = calculateMST(clusterCenters);

                // Update the graph panel with a new visualizer
                visualizer = new MetroClusterVisualizer(planner.getData(), clusters, clusterCenters, finalStations, mstEdges);
                graphPanel.removeAll();
                graphPanel.add(visualizer, BorderLayout.CENTER);
                graphPanel.revalidate();
                graphPanel.repaint();

                centralityButton.setEnabled(true);
                JOptionPane.showMessageDialog(frame, "Clustering Applied and Graph Displayed Successfully");
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(frame, "Error in clustering: " + ex.getMessage());
            }
        });

        // Action listener for betweenness centrality
        centralityButton.addActionListener(e -> {
            try {
                finalStations = calculateFinalStations(clusterCenters);
                System.out.println(finalStations);
                visualizer.updateFinalStations(finalStations);
                graphPanel.repaint();

                JOptionPane.showMessageDialog(frame,
                        "Final Metro Stations Selected\nNumber of stations: " + finalStations.size());
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(frame, "Error in centrality calculation: " + ex.getMessage());
            }
        });

        // Add components to frame
        frame.add(fileButton);
        frame.add(preprocessButton);
        frame.add(clusterButton);
        frame.add(centralityButton);
        frame.add(graphPanel);

        frame.setVisible(true);
    }

    private Map<Integer, Point2D> calculateClusterCenters(List<Map<String, Object>> data, int[] clusters) {
        Map<Integer, List<Point2D>> clusterPoints = new HashMap<>();
        Map<Integer, Point2D> centers = new HashMap<>();

        // Group points by cluster
        for (int i = 0; i < clusters.length; i++) {
            int cluster = clusters[i];
            if (cluster != -1) {
                double lon = Double.parseDouble(data.get(i).get("longitude").toString());
                double lat = Double.parseDouble(data.get(i).get("latitude").toString());
                clusterPoints.computeIfAbsent(cluster, k -> new ArrayList<>())
                        .add(new Point2D.Double(lon, lat));
            }
        }

        // Calculate center for each cluster
        for (Map.Entry<Integer, List<Point2D>> entry : clusterPoints.entrySet()) {
            List<Point2D> points = entry.getValue();
            double avgX = points.stream().mapToDouble(p -> p.getX()).average().orElse(0.0);
            double avgY = points.stream().mapToDouble(p -> p.getY()).average().orElse(0.0);
            centers.put(entry.getKey(), new Point2D.Double(avgX, avgY));
        }

        return centers;
    }

    private Set<Integer> calculateFinalStations(Map<Integer, Point2D> clusterCenters) {
        Set<Integer> selectedStations = new HashSet<>();

        // Create a graph representation
        Map<Integer, Map<Integer, Double>> graph = new HashMap<>();
        for (Integer node1 : clusterCenters.keySet()) {
            graph.put(node1, new HashMap<>());
            for (Integer node2 : clusterCenters.keySet()) {
                if (!node1.equals(node2)) {
                    Point2D p1 = clusterCenters.get(node1);
                    Point2D p2 = clusterCenters.get(node2);
                    double distance = p1.distance(p2);
                    graph.get(node1).put(node2, distance);
                }
            }
        }

        // Calculate betweenness centrality
        Map<Integer, Double> betweenness = calculateBetweennessCentrality(graph);

        // Select top stations based on betweenness centrality
        List<Map.Entry<Integer, Double>> sortedStations = new ArrayList<>(betweenness.entrySet());
        sortedStations.sort((a, b) -> b.getValue().compareTo(a.getValue()));

        // Select top 30% of stations
        int numStations = (int) Math.ceil(clusterCenters.size() * 0.3);
        for (int i = 0; i < numStations && i < sortedStations.size(); i++) {
            selectedStations.add(sortedStations.get(i).getKey());
        }

        return selectedStations;
    }

    private Map<Integer, Double> calculateBetweennessCentrality(Map<Integer, Map<Integer, Double>> graph) {
        Map<Integer, Double> betweenness = new HashMap<>();
        for (Integer v : graph.keySet()) {
            betweenness.put(v, 0.0);
        }

        // For each pair of nodes
        for (Integer s : graph.keySet()) {
            for (Integer t : graph.keySet()) {
                if (s.equals(t)) continue;

                // Calculate shortest paths using Dijkstra's algorithm
                Map<Integer, Double> distances = new HashMap<>();
                Map<Integer, List<Integer>> predecessors = new HashMap<>();
                floydwarshallShortestPaths(graph, s, distances, predecessors);

                // Update betweenness values
                updateBetweenness(s, t, predecessors, betweenness);
            }
        }

        return betweenness;
    }

    private void floydwarshallShortestPaths(
            Map<Integer, Map<Integer, Double>> graph,
            Integer source,
            Map<Integer, Double> distances,
            Map<Integer, List<Integer>> predecessors) {

        // Get all unique vertices
        Set<Integer> vertices = new HashSet<>(graph.keySet());

        // Initialize distance matrix
        Map<Integer, Map<Integer, Double>> dist = new HashMap<>();
        Map<Integer, Map<Integer, Integer>> next = new HashMap<>();

        // Initialize distances and predecessors
        for (Integer v : vertices) {
            dist.put(v, new HashMap<>());
            next.put(v, new HashMap<>());
            predecessors.put(v, new ArrayList<>());

            for (Integer w : vertices) {
                if (v.equals(w)) {
                    dist.get(v).put(w, 0.0);
                } else if (graph.containsKey(v) && graph.get(v).containsKey(w)) {
                    dist.get(v).put(w, graph.get(v).get(w));
                    next.get(v).put(w, w);
                } else {
                    dist.get(v).put(w, Double.MAX_VALUE);
                }
            }
        }

        // Floyd-Warshall algorithm
        for (Integer k : vertices) {
            for (Integer i : vertices) {
                for (Integer j : vertices) {
                    // Prevent potential overflow
                    if (dist.get(i).get(k) != Double.MAX_VALUE &&
                            dist.get(k).get(j) != Double.MAX_VALUE) {
                        double newDistance = dist.get(i).get(k) + dist.get(k).get(j);

                        if (newDistance < dist.get(i).get(j)) {
                            dist.get(i).put(j, newDistance);
                            next.get(i).put(j, next.get(i).get(k));
                        }
                    }
                }
            }
        }

        // Populate distances and predecessors for the source
        for (Integer v : vertices) {
            if (!v.equals(source)) {
                double pathDistance = dist.get(source).get(v);
                distances.put(v, pathDistance);

                // Reconstruct path
                if (pathDistance != Double.MAX_VALUE) {
                    Integer current = source;
                    while (!current.equals(v)) {
                        current = next.get(current).get(v);
                        predecessors.get(v).add(0, current);
                    }
                }
            }
        }
    }

    private void updateBetweenness(
            Integer source,
            Integer target,
            Map<Integer, List<Integer>> predecessors,
            Map<Integer, Double> betweenness) {

        Set<Integer> visited = new HashSet<>();
        Stack<Integer> stack = new Stack<>();
        stack.push(target);

        while (!stack.isEmpty()) {
            Integer v = stack.pop();
            if (!v.equals(source)) {
                for (Integer pred : predecessors.get(v)) {
                    if (!visited.contains(pred)) {
                        stack.push(pred);
                        betweenness.put(pred, betweenness.get(pred) + 1.0);
                    }
                }
            }
            visited.add(v);
        }
    }

    private Set<Line2D> calculateMST(Map<Integer, Point2D> clusterCenters) {
        Set<Line2D> mstEdges = new HashSet<>();
        Map<Integer, Boolean> visited = new HashMap<>();
        PriorityQueue<Edge> pq = new PriorityQueue<>(Comparator.comparingDouble(e -> e.weight));

        // Only consider final stations
        Map<Integer, Point2D> finalStationCenters = new HashMap<>();
        for (Integer id : finalStations) {
            if (clusterCenters.containsKey(id)) {
                finalStationCenters.put(id, clusterCenters.get(id));
            }
        }

        // If there are no final stations or only one, return empty set
        if (finalStationCenters.size() < 2) {
            return mstEdges;
        }

        // Initialize visited map
        for (Integer node : finalStationCenters.keySet()) {
            visited.put(node, false);
        }

        // Start from an arbitrary node
        Integer startNode = finalStationCenters.keySet().iterator().next();
        visited.put(startNode, true);

        // Add all edges from the start node to the priority queue
        for (Integer neighbor : finalStationCenters.keySet()) {
            if (!neighbor.equals(startNode)) {
                double distance = finalStationCenters.get(startNode).distance(finalStationCenters.get(neighbor));
                pq.add(new Edge(startNode, neighbor, distance));
            }
        }

        // Prim's algorithm
        while (!pq.isEmpty()) {
            Edge edge = pq.poll();
            if (!visited.get(edge.to)) {
                visited.put(edge.to, true);
                mstEdges.add(new Line2D.Double(finalStationCenters.get(edge.from),
                        finalStationCenters.get(edge.to)));

                // Add new edges to the priority queue
                for (Integer neighbor : finalStationCenters.keySet()) {
                    if (!visited.get(neighbor)) {
                        double distance = finalStationCenters.get(edge.to).distance(finalStationCenters.get(neighbor));
                        pq.add(new Edge(edge.to, neighbor, distance));
                    }
                }
            }
        }

        return mstEdges;
    }

    private static class Edge {
        int from, to;
        double weight;

        Edge(int from, int to, double weight) {
            this.from = from;
            this.to = to;
            this.weight = weight;
        }
    }

    public static void main(String[] args) {
        new MetroStationPlannerGUI();
    }
}

class MetroClusterVisualizer extends JPanel implements MouseWheelListener, MouseMotionListener {
    private List<Map<String, Object>> data;
    private int[] clusters;
    private Map<Integer, Color> clusterColors;
    private Map<Integer, Point2D> clusterCenters;
    private Set<Integer> finalStations;
    private Set<Line2D> mstEdges;

    private double zoomFactor = 1.0;
    private int translateX = 0, translateY = 0;
    private int prevX, prevY;

    public MetroClusterVisualizer(List<Map<String, Object>> data, int[] clusters,
                                  Map<Integer, Point2D> clusterCenters, Set<Integer> finalStations, Set<Line2D> mstEdges) {
        this.data = data;
        this.clusters = clusters;
        this.clusterCenters = clusterCenters;
        this.finalStations = finalStations;
        this.clusterColors = assignClusterColors();
        this.mstEdges = mstEdges;

        addMouseWheelListener(this);
        addMouseMotionListener(this);
    }

    public void updateFinalStations(Set<Integer> finalStations) {
        this.finalStations = finalStations;

        // Create MST for final stations only
        Set<Line2D> newMstEdges = new HashSet<>();
        if (finalStations.size() < 2) {
            this.mstEdges = newMstEdges;
            repaint();
            return;
        }

        // Use Prim's algorithm for MST
        Map<Integer, Boolean> visited = new HashMap<>();
        for (Integer station : finalStations) {
            visited.put(station, false);
        }

        // Start with first station
        Integer start = finalStations.iterator().next();
        visited.put(start, true);

        // Create edges until all stations are connected
        while (visited.containsValue(false)) {
            double minDistance = Double.MAX_VALUE;
            Line2D bestEdge = null;
            Integer nextStation = null;

            // Find shortest edge to unvisited station
            for (Integer from : finalStations) {
                if (!visited.get(from)) continue;

                for (Integer to : finalStations) {
                    if (visited.get(to)) continue;

                    Point2D fromPoint = clusterCenters.get(from);
                    Point2D toPoint = clusterCenters.get(to);
                    double distance = fromPoint.distance(toPoint);

                    if (distance < minDistance) {
                        minDistance = distance;
                        bestEdge = new Line2D.Double(fromPoint, toPoint);
                        nextStation = to;
                    }
                }
            }

            if (bestEdge != null) {
                newMstEdges.add(bestEdge);
                visited.put(nextStation, true);
            }
        }

        this.mstEdges = newMstEdges;
        System.out.println("Updated final stations: " + finalStations);
        System.out.println("Number of MST paths: " + mstEdges.size());
        repaint();
    }


    private Map<Integer, Color> assignClusterColors() {
        Map<Integer, Color> colorMap = new HashMap<>();
        Random rand = new Random(42); // Fixed seed for consistent colors
        for (int cluster : clusters) {
            if (cluster == -1) continue;
            if (!colorMap.containsKey(cluster)) {
                colorMap.put(cluster, new Color(rand.nextInt(255), rand.nextInt(255), rand.nextInt(255)));
            }
        }
        return colorMap;
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        int legendWidth = 150;
        int graphWidth = getWidth() - legendWidth - 20;
        int graphHeight = getHeight() - 100;

        // Find bounds
        double minLon = Double.MAX_VALUE, maxLon = Double.MIN_VALUE;
        double minLat = Double.MAX_VALUE, maxLat = Double.MIN_VALUE;
        for (Map<String, Object> point : data) {
            double lon = Double.parseDouble(point.get("longitude").toString());
            double lat = Double.parseDouble(point.get("latitude").toString());
            minLon = Math.min(minLon, lon);
            maxLon = Math.max(maxLon, lon);
            minLat = Math.min(minLat, lat);
            maxLat = Math.max(maxLat, lat);
        }

        // Transform setup
        AffineTransform transform = new AffineTransform();
        transform.translate(translateX, translateY);
        transform.scale(zoomFactor, zoomFactor);
        g2.setTransform(transform);

        // Draw axes
        g2.setColor(Color.BLACK);
        g2.draw(new Line2D.Double(50, graphHeight + 50, graphWidth + 50, graphHeight + 50));
        g2.draw(new Line2D.Double(50, 50, 50, graphHeight + 50));
        g2.drawString("Longitude", graphWidth / 2 + 50, graphHeight + 70);
        g2.drawString("Latitude", 10, graphHeight / 2 + 50);

        // Draw cluster points and centers
        for (int i = 0; i < data.size(); i++) {
            int cluster = clusters[i];
            if (cluster != -1) {
                double lon = Double.parseDouble(data.get(i).get("longitude").toString());
                double lat = Double.parseDouble(data.get(i).get("latitude").toString());

                int x = (int) ((lon - minLon) / (maxLon - minLon) * graphWidth) + 50;
                int y = graphHeight + 50 - (int) ((lat - minLat) / (maxLat - minLat) * graphHeight);

                g2.setColor(clusterColors.get(cluster));
                g2.fill(new Ellipse2D.Double(x - 3, y - 3, 6, 6));
            }
        }

        // Draw cluster centers and final stations
        for (Map.Entry<Integer, Point2D> entry : clusterCenters.entrySet()) {
            int clusterId = entry.getKey();
            Point2D center = entry.getValue();

            int x = (int) ((center.getX() - minLon) / (maxLon - minLon) * graphWidth) + 50;
            int y = graphHeight + 50 - (int) ((center.getY() - minLat) / (maxLat - minLat) * graphHeight);

            // Draw larger circles for cluster centers
            g2.setColor(Color.BLACK);
            g2.setStroke(new BasicStroke(2.0f));
            g2.draw(new Ellipse2D.Double(x - 8, y - 8, 16, 16));

            // Highlight final stations with a star shape
            if (finalStations.contains(clusterId)) {
                drawStar(g2, x, y, 15);
                g2.setColor(Color.RED);
                g2.setFont(new Font("Arial", Font.BOLD, 12));
                g2.drawString("Station " + clusterId, x + 20, y);
            }
        }

        // Draw paths between final stations
        drawPathsBetweenFinalStations(g2, graphWidth, graphHeight, minLon, maxLon, minLat, maxLat);

        // Draw legend
        drawLegend(g2, graphWidth + 70, 50);
    }

    private void drawPathsBetweenFinalStations(Graphics2D g2, int graphWidth, int graphHeight,
                                               double minLon, double maxLon, double minLat, double maxLat) {
        // Set prominent style for the paths
        g2.setColor(new Color(0, 100, 200)); // Darker blue
        g2.setStroke(new BasicStroke(2.5f)); // Thicker lines

        for (Line2D edge : mstEdges) {
            Point2D p1 = edge.getP1();
            Point2D p2 = edge.getP2();

            int x1 = (int) ((p1.getX() - minLon) / (maxLon - minLon) * graphWidth) + 50;
            int y1 = graphHeight + 50 - (int) ((p1.getY() - minLat) / (maxLat - minLat) * graphHeight);
            int x2 = (int) ((p2.getX() - minLon) / (maxLon - minLon) * graphWidth) + 50;
            int y2 = graphHeight + 50 - (int) ((p2.getY() - minLat) / (maxLat - minLat) * graphHeight);

            g2.draw(new Line2D.Double(x1, y1, x2, y2));
        }
    }

    private void drawStar(Graphics2D g2, int x, int y, int size) {
        int numPoints = 5;
        int[] xPoints = new int[numPoints * 2];
        int[] yPoints = new int[numPoints * 2];

        double angle = Math.PI / numPoints;

        for (int i = 0; i < numPoints * 2; i++) {
            double r = (i % 2 == 0) ? size : size/2;
            xPoints[i] = (int) (x + Math.cos(i * angle) * r);
            yPoints[i] = (int) (y + Math.sin(i * angle) * r);
        }

        g2.setColor(Color.RED);
        g2.setStroke(new BasicStroke(2.0f));
        g2.drawPolygon(xPoints, yPoints, numPoints * 2);
        g2.setColor(new Color(255, 200, 200));
        g2.fillPolygon(xPoints, yPoints, numPoints * 2);
    }

    private void drawLegend(Graphics2D g2, int startX, int startY) {
        g2.setColor(Color.BLACK);
        g2.setFont(new Font("Arial", Font.BOLD, 14));
        g2.drawString("Legend:", startX, startY);

        int yOffset = 30;

        // Regular points
        g2.setColor(Color.BLACK);
        g2.fillOval(startX, startY + yOffset, 6, 6);
        g2.drawString("Data Points", startX + 20, startY + yOffset + 6);

        // Cluster centers
        yOffset += 25;
        g2.setColor(Color.BLACK);
        g2.setStroke(new BasicStroke(2.0f));
        g2.drawOval(startX - 2, startY + yOffset - 2, 12, 12);
        g2.drawString("Cluster Centers", startX + 20, startY + yOffset + 6);

        // Final stations
        yOffset += 25;
        drawStar(g2, startX + 5, startY + yOffset + 5, 10);
        g2.setColor(Color.BLACK);
        g2.drawString("Metro Stations", startX + 20, startY + yOffset + 6);

        // Cluster colors
        yOffset += 30;
        g2.drawString("Clusters:", startX, startY + yOffset);
        for (Map.Entry<Integer, Color> entry : clusterColors.entrySet()) {
            yOffset += 20;
            g2.setColor(entry.getValue());
            g2.fillRect(startX, startY + yOffset, 15, 15);
            g2.setColor(Color.BLACK);
            g2.drawString("Cluster " + entry.getKey(), startX + 20, startY + yOffset + 12);
        }
    }

    @Override
    public void mouseWheelMoved(MouseWheelEvent e) {
        if (e.getWheelRotation() < 0) {
            zoomFactor *= 1.1;
        } else {
            zoomFactor /= 1.1;
        }
        repaint();
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        int deltaX = e.getX() - prevX;
        int deltaY = e.getY() - prevY;

        translateX += deltaX;
        translateY += deltaY;

        prevX = e.getX();
        prevY = e.getY();

        repaint();
    }

    @Override
    public void mouseMoved(MouseEvent e) {
        prevX = e.getX();
        prevY = e.getY();
    }
}