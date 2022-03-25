package org.locationtech.jts.io.twkb;

import org.locationtech.jts.geom.CoordinateSequence;
import org.locationtech.jts.geom.CoordinateSequenceFilter;

class BoundsExtractor implements CoordinateSequenceFilter {

    private final int dimensions;

    double[] ordinates = new double[]{ //
        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, // note, Double.MIN_VALUE is positive
        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, //
        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, //
        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY//
    };

    BoundsExtractor(int dimensions) {
        this.dimensions = dimensions;
    }

    public @Override
    void filter(final CoordinateSequence seq, final int coordIndex) {
        for (int ordinateIndex = 0; ordinateIndex < dimensions; ordinateIndex++) {
            final double ordinate = seq.getOrdinate(coordIndex, ordinateIndex);
            final int minIndex = 2 * ordinateIndex;
            final int maxIndex = minIndex + 1;
            ordinates[minIndex] = Math.min(ordinates[minIndex], ordinate);
            ordinates[maxIndex] = Math.max(ordinates[maxIndex], ordinate);
        }
    }

    @Override
    public boolean isDone() {
        return false;
    }

    @Override
    public boolean isGeometryChanged() {
        return false;
    }
}
