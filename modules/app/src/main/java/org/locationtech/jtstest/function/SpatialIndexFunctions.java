/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

package org.locationtech.jtstest.function;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.GeometryFilter;
import org.locationtech.jts.index.SpatialIndex;
import org.locationtech.jts.index.hprtree.HPRtree;
import org.locationtech.jts.index.kdtree.KdTree;
import org.locationtech.jts.index.quadtree.Quadtree;
import org.locationtech.jts.index.strtree.AbstractNode;
import org.locationtech.jts.index.strtree.Boundable;
import org.locationtech.jts.index.strtree.GeometryItemDistance;
import org.locationtech.jts.index.strtree.STRtree;


public class SpatialIndexFunctions
{
  public static Geometry kdTreeQuery(Geometry pts, Geometry queryEnv, double tolerance)
  {
    KdTree index = buildKdTree(pts, tolerance);
    // if no query env provided query everything inserted 
    if (queryEnv == null) queryEnv = pts;
    List result = index.query(queryEnv.getEnvelopeInternal());
    Coordinate[] resultCoords = KdTree.toCoordinates(result);
    return pts.getFactory().createMultiPoint(resultCoords);
  }

  public static Geometry kdTreeQueryRepeated(Geometry pts, Geometry queryEnv, double tolerance)
  {
    KdTree index = buildKdTree(pts, tolerance);
    // if no query env provided query everything inserted 
    if (queryEnv == null) queryEnv = pts;
    List result = index.query(queryEnv.getEnvelopeInternal());
    Coordinate[] resultCoords = KdTree.toCoordinates(result, true);
    return pts.getFactory().createMultiPoint(resultCoords);
  }

  private static KdTree buildKdTree(Geometry geom, double tolerance) {
    final KdTree index = new KdTree(tolerance);
    Coordinate[] pt = geom.getCoordinates();
    for (int i = 0; i < pt.length; i++) {
      index.insert(pt[i]);
    }
    return index;
  }
  
  public static Geometry strTreeBounds(Geometry geoms)
  {
    STRtree index = new STRtree();
    loadIndex(geoms, index);
    List bounds = new ArrayList();
    addBounds(index.getRoot(), bounds, geoms.getFactory());
    return geoms.getFactory().buildGeometry(bounds);
  }

  private static void addBounds(Boundable bnd, List bounds,
      GeometryFactory factory) {
    // don't include bounds of leaf nodes
    if (! (bnd instanceof AbstractNode)) return;
    
    Envelope env = (Envelope) bnd.getBounds();
    bounds.add(factory.toGeometry(env));
    if (bnd instanceof AbstractNode) {
      AbstractNode node = (AbstractNode) bnd;
      List children = node.getChildBoundables();
      for (Iterator i = children.iterator(); i.hasNext(); ) {
        Boundable child = (Boundable) i.next();
        addBounds(child, bounds, factory);
      }
    }
  }

  public static Geometry hprTreeQuery(Geometry geoms, Geometry queryEnv)
  {
    HPRtree index = new HPRtree();
    loadIndex(geoms, index);
    // if no query env provided query everything inserted 
    if (queryEnv == null) queryEnv = geoms;
    List result = index.query(queryEnv.getEnvelopeInternal());
    return geoms.getFactory().buildGeometry(result);
  }

  private static void loadIndex(Geometry geom, SpatialIndex index) {
    geom.apply(new GeometryFilter() {

      public void filter(Geometry geom) {
        // only insert atomic geometries
        if (geom instanceof GeometryCollection) return;
        index.insert(geom.getEnvelopeInternal(), geom);
      }
      
    });
  }

  public static Geometry hprTreeBounds(Geometry geoms)
  {
    HPRtree index = new HPRtree();
    loadIndex(geoms, index);
    index.build();
    Envelope[] bounds = index.getBounds();
    Geometry[] polys = new Geometry[bounds.length];
    int i = 0;
    for (Envelope env : bounds) {
      polys[i++] = geoms.getFactory().toGeometry(env);
    }
    return geoms.getFactory().createGeometryCollection(polys);
  }
  
  public static Geometry strTreeQuery(Geometry geoms, Geometry queryEnv)
  {
    STRtree index = new STRtree();
    loadIndex(geoms, index);
    // if no query env provided query everything inserted 
    if (queryEnv == null) queryEnv = geoms;
    List result = index.query(queryEnv.getEnvelopeInternal());
    return geoms.getFactory().buildGeometry(result);
  }
  
  public static Geometry strTreeNN(Geometry geoms, Geometry geom)
  {
    STRtree index = new STRtree();
    loadIndex(geoms, index);
    Object result = index.nearestNeighbour(geom.getEnvelopeInternal(), geom, new GeometryItemDistance());
    return (Geometry) result;
  }

  public static Geometry strTreeNNInSet(Geometry geoms)
  {
    STRtree index = new STRtree();
    loadIndex(geoms, index);
    Object[] result = index.nearestNeighbour(new GeometryItemDistance());
    Geometry[] resultGeoms = new Geometry[] { (Geometry) result[0], (Geometry) result[1] };
    return geoms.getFactory().createGeometryCollection(resultGeoms);
  }

  public static Geometry strTreeNNk(Geometry geoms, Geometry geom, int k)
  {
    STRtree index = new STRtree();
    loadIndex(geoms, index);
    Object[] knnObjects = index.nearestNeighbour(geom.getEnvelopeInternal(), geom, new GeometryItemDistance(), k);
    List knnGeoms = new ArrayList(Arrays.asList(knnObjects));
    Geometry geometryCollection = geoms.getFactory().buildGeometry(knnGeoms);
    return geometryCollection;
  }
  
  public static Geometry quadTreeQuery(Geometry geoms, Geometry queryEnv)
  {
    Quadtree index = buildQuadtree(geoms);
    // if no query env provided query everything inserted 
    if (queryEnv == null) queryEnv = geoms;
    List result = index.query(queryEnv.getEnvelopeInternal());
    return geoms.getFactory().buildGeometry(result);
  }

  private static Quadtree buildQuadtree(Geometry geom) {
    final Quadtree index = new Quadtree();
    geom.apply(new GeometryFilter() {

      public void filter(Geometry geom) {
        // only insert atomic geometries
        if (geom instanceof GeometryCollection) return;
        index.insert(geom.getEnvelopeInternal(), geom);
      }
      
    });
    return index;
  }
}
