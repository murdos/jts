/*
 * Copyright (c) 2018 James Hughes, 2019 Gabriel Roldan
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */
package org.locationtech.jts.io.twkb;

import static org.locationtech.jts.io.twkb.TWKBHeader.GeometryType.POINT;

import java.io.ByteArrayOutputStream;
import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Objects;

import org.locationtech.jts.geom.CoordinateSequence;
import org.locationtech.jts.geom.CoordinateSequenceFilter;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.io.twkb.TWKBHeader.GeometryType;

/**
 * <pre>
 * {@code
 * 
 * twkb                    := <header> <geometry_body>
 * header                  := <type_and_precision> <metadata_header> [extended_dimensions_header] [geometry_body_size]
 * type_and_precision      := byte := <type_mask OR precision_mask>)
 * type_mask               := <ubyte> (0b0000XXXX -> 1=point, 2=linestring, 3=polygon, 4=multipoint, 
 *                                     5=multilinestring, 6=multipolygon, 7=geometry collection)
 * precision_mask          := <signed byte> (zig-zag encoded 4-bit signed integer, 0bXXXX0000. Number of base-10 decimal places 
 *                                           stored. A positive retaining information to the right of the decimal place, negative 
 *                                           rounding up to the left of the decimal place)  
 * metadata_header := byte := <bbox_flag OR  size_flag OR idlist_flag OR extended_precision_flag OR empty_geometry_flag>
 * bbox_flag               := 0b00000001
 * size_flag               := 0b00000010
 * idlist_flag             := 0b00000100
 * extended_precision_flag := 0b00001000
 * empty_geometry_flag     := 0b00010000
 * 
 * # extended_dimensions_header present iif extended_precision_flag == 1
 * extended_dimensions_header  := byte := <Z_coordinates_presence_flag OR M_coordinates_presence_flag OR Z_precision OR M_precision>
 * Z_coordinates_presence_flag := 0b00000001 
 * M_coordinates_presence_flag := 0b00000010
 * Z_precision                 := 0b000XXX00 3-bit unsigned integer using bits 3-5 
 * M_precision                 := 0bXXX00000 3-bit unsigned integer using bits 6-8
 * 
 * # geometry_body_size present iif size_flag == 1 
 * geometry_body_size := uint32 # size in bytes of <geometry_body>
 * 
 * # geometry_body present iif empty_geometry_flag == 0
 * geometry_body := [bounds] [idlist] <geometry>
 * # bounds present iff bbox_flag == 1 
 * # 2 signed varints per dimension. i.e.:
 * # [xmin, deltax, ymin, deltay]                              iif Z_coordinates_presence_flag == 0 AND M_coordinates_presence_flag == 0
 * # [xmin, deltax, ymin, deltay, zmin, deltaz]                iif Z_coordinates_presence_flag == 1 AND M_coordinates_presence_flag == 0
 * # [xmin, deltax, ymin, deltay, zmin, deltaz, mmin, deltam]  iif Z_coordinates_presence_flag == 1 AND M_coordinates_presence_flag == 1
 * # [xmin, deltax, ymin, deltay, mmin, deltam]                iif Z_coordinates_presence_flag == 0 AND M_coordinates_presence_flag == 1
 * bounds          := sint32[4] | sint32[6] | sint32[8] 
 * geometry        := point | linestring | polygon | multipoint | multilinestring | multipolygon | geomcollection
 * point           := sint32[dimension]
 * linestring      := <npoints:uint32> [point[npoints]]
 * polygon         := <nrings:uint32> [linestring]
 * multipoint      := <nmembers:uint32> [idlist:<sint32[nmembers]>] [point[nmembers]]
 * multilinestring := <nmembers:uint32> [idlist:<sint32[nmembers]>] [linestring[nmembers]]
 * multipolygon    := <nmembers:uint32> [idlist:<sint32[nmembers]>] [polygon[nmembers]]
 * geomcollection  := <nmembers:uint32> [idlist:<sint32[nmembers]>] [twkb[nmembers]]
 * 
 * uint32 := <Unsigned variable-length encoded integer>
 * sint32 := <Signed variable-length, zig-zag encoded integer>
 * byte := <Single octect>
 * 
 * }
 * </pre>
 */
public class TWKBWriter {

    private TWKBHeader paramsHeader = new TWKBHeader()
        .setXyPrecision(7)
        .setZPrecision(0)
        .setMPrecision(0);

    private boolean optimizedEncoding = true;

    /**
     * Number of base-10 decimal places stored for X and Y dimensions.
     * <p>
     * A positive retaining information to the right of the decimal place, negative rounding up to
     * the left of the decimal place).
     * <p>
     * Defaults to {@code 7}
     */
    public TWKBWriter setXYPrecision(int xyprecision) {
        if (xyprecision < -7 || xyprecision > 7) {
            throw new IllegalArgumentException(
                    "X/Z precision cannot be greater than 7 or less than -7");
        }
        paramsHeader = paramsHeader.setXyPrecision(xyprecision);
        return this;
    }

    public TWKBWriter setEncodeZ(boolean includeZDimension) {
        paramsHeader = paramsHeader.setHasZ(includeZDimension);
        return this;
    }

    public TWKBWriter setEncodeM(boolean includeMDimension) {
        paramsHeader = paramsHeader.setHasM(includeMDimension);
        return this;
    }

    /**
     * Number of base-10 decimal places stored for Z dimension.
     * <p>
     * A positive retaining information to the right of the decimal place, negative rounding up to
     * the left of the decimal place).
     * <p>
     * Defaults to {@code 0}
     */
    public TWKBWriter setZPrecision(int zprecision) {
        if (zprecision < 0 || zprecision > 7) {
            throw new IllegalArgumentException("Z precision cannot be negative or greater than 7");
        }
        paramsHeader = paramsHeader.setZPrecision(zprecision);
        return this;
    }

    /**
     * Number of base-10 decimal places stored for M dimension.
     * <p>
     * A positive retaining information to the right of the decimal place, negative rounding up to
     * the left of the decimal place).
     * <p>
     * Defaults to {@code 0}
     */
    public TWKBWriter setMPrecision(int mprecision) {
        if (mprecision < 0 || mprecision > 7) {
            throw new IllegalArgumentException("M precision cannot be negative or greater than 7");
        }
        paramsHeader = paramsHeader.setMPrecision(mprecision);
        return this;
    }

    /**
     * Whether the generated TWKB should include the size in bytes of the geometry.
     */
    public TWKBWriter setIncludeSize(boolean includeSize) {
        paramsHeader = paramsHeader.setHasSize(includeSize);
        return this;
    }

    /**
     * Whether the generated TWKB should include a Bounding Box for the geometry.
     */
    public TWKBWriter setIncludeBbox(boolean includeBbox) {
        paramsHeader = paramsHeader.setHasBBOX(includeBbox);
        return this;
    }

    /**
     * Enables or disables the following optimizations at encoding time, defaults to {@code true}:
     * <ul>
     * <li>The {@code xy}, {@code m}, and {@code z} precision of an {@link Geometry#isEmpty() empty}
     * geometry are set to {@code 0}, despite the values of {@link #setXYPrecision(int)},
     * {@link #setZPrecision(int)}, and {@link #setMPrecision(int)}
     * <li>BBOX is not encoded for {@link Point} geometries, despite {@link #setIncludeBbox(boolean)} being
     * {@code true}, and the {@code hasBBOX} flag is encoded as {@code false}
     * <li>{@link TWKBHeader#hasZ()} is forcedly encoded as {@code false} if the input geometry has no Z
     * dimension (as per {@link CoordinateSequence#hasZ()}), even if {@code hasZ() == true}, and
     * consequently {@link TWKBHeader#zPrecision()} is not set (defaulting to 0 if there's an M dimension).
     * This fixes a possible mismatch in the metadata where a dimension precision may be present in
     * the header where the geometry doesn't really have that dimension.
     * <li>{@link TWKBHeader#hasM()} is forcedly encoded as {@code false} if the input geometry has no M
     * dimension (as per {@link CoordinateSequence#hasM()}), even if {@code hasM() == true}, and
     * consequently {@link TWKBHeader#mPrecision()} is not set (defaulting to 0 if there's a Z dimension).
     * This fixes a possible mismatch in the metadata where a dimension precision may be present in
     * the header where the geometry doesn't really have that dimension.
     * <li>For a <strong>Polygon</strong> geometry, all its {@link LinearRing}s get their last
     * coordinate removed, as long as that wouldn't result in an invalid geometry. Following the
     * specification's suggestion, the last coordinate of a {@link LinearRing} is removed and
     * therefore it's encoded as a {@link LineString} with one less coordinate. When parsing, such
     * line strings are automatically closed to form a {@code LinearRing}
     * <li>Consecutive coordinates of a {@link CoordinateSequence} that, due to loss of precision,
     * collapse to a single point, are merged into a single coordinate, as long as that doesn't
     * result in an invalid geometry. For example,
     * {@code LINESTRING(0.1 0.1, 0.15 0.15, 1.2 1.2, 1.29 1.29)}, encoded with an {@code xy}
     * precision of {@code 1} decimal places, is collapsed to {@code LINESTRING(0.1 0.1, 1.2 1.2)}
     * instead of {@code LINESTRING(0.1 0.1, 0.1 0.1, 1.2 1.2, 1.2 1.2)}
     * </ul>
     * Disabling these optimizations make the encoding consistent with PostGIS TWKB encoding.
     */
    public TWKBWriter setOptimizedEncoding(boolean optimizedEncoding) {
        this.optimizedEncoding = optimizedEncoding;
        return this;
    }

    public byte[] write(Geometry geom) {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        try {
            write(geom, out);
        } catch (IOException ex) {
            throw new RuntimeException("Unexpected IOException caught: " + ex.getMessage(), ex);
        }
        return out.toByteArray();
    }

    public void write(Geometry geom, OutputStream out) throws IOException {
        write(geom, (DataOutput) new DataOutputStream(out));
    }

    public void write(Geometry geom, DataOutput out) throws IOException {
        writeInternal(geom, out);
    }

    final /* @VisibleForTesting */ TWKBHeader writeInternal(Geometry geom, DataOutput out)
            throws IOException {
        Objects.requireNonNull(geom, "geometry is null");
        return write(geom, TWKBOutputStream.of(out), paramsHeader);
    }

    private TWKBHeader write(Geometry geometry, TWKBOutputStream out, TWKBHeader params)
        throws IOException {
        return write(geometry, out, params, false);
    }

    private TWKBHeader write(Geometry geometry, TWKBOutputStream out, TWKBHeader params,
        boolean forcePreserveHeaderDimensions) throws IOException {
        Objects.requireNonNull(geometry, "Geometry is null");
        Objects.requireNonNull(out, "DataOutput is null");
        Objects.requireNonNull(params, "TWKBHeader is null");

        TWKBHeader header = prepareHeader(geometry, new TWKBHeader(params), forcePreserveHeaderDimensions);

        if (header.hasSize()) {
            BufferedTKWBOutputStream bufferedBody = BufferedTKWBOutputStream.create();
            writeGeometryBody(geometry, bufferedBody, header);
            int bodySize = bufferedBody.size();
            header = header.setGeometryBodySize(bodySize);
            writeHeaderTo(header, out);
            bufferedBody.writeTo(out);
        } else {
            writeHeaderTo(header, out);
            writeGeometryBody(geometry, out, header);
        }
        return header;
    }

    static void writeHeaderTo(TWKBHeader header, DataOutput out) throws IOException {
        writeHeaderTo(header, TWKBOutputStream.of(out));
    }

    private static void writeHeaderTo(TWKBHeader header, TWKBOutputStream out) throws IOException {
        Objects.requireNonNull(out);
        final int typeAndPrecisionHeader;
        final int metadataHeader;
        {
            final int geometryType = header.geometryType().getValue();
            final int precisionHeader = Varint.zigZagEncode(header.xyPrecision()) << 4;
            typeAndPrecisionHeader = precisionHeader | geometryType;

            metadataHeader = (header.hasBBOX() ? 0b00000001 : 0) //
                | (header.hasSize() ? 0b00000010 : 0)//
                | (header.hasIdList() ? 0b00000100 : 0)//
                | (header.hasExtendedPrecision() ? 0b00001000 : 0)//
                | (header.isEmpty() ? 0b00010000 : 0);
        }
        out.writeByte(typeAndPrecisionHeader);
        out.writeByte(metadataHeader);
        if (header.hasExtendedPrecision()) {
            // final int extendedDimsHeader = in.readByte() & 0xFF;
            // hasZ = (extendedDimsHeader & 0b00000001) > 0;
            // hasM = (extendedDimsHeader & 0b00000010) > 0;
            // zprecision = (extendedDimsHeader & 0b00011100) >> 2;
            // mprecision = (extendedDimsHeader & 0b11100000) >> 5;

            int extendedDimsHeader = (header.hasZ() ? 0b00000001 : 0) | (header.hasM() ? 0b00000010 : 0);
            extendedDimsHeader |= header.zPrecision() << 2;
            extendedDimsHeader |= header.mPrecision() << 5;

            out.writeByte(extendedDimsHeader);
        }
        if (header.hasSize()) {
            out.writeUnsignedVarInt(header.geometryBodySize());
        }
    }

    private TWKBHeader prepareHeader(Geometry geometry, TWKBHeader params,
        boolean forcePreserveHeaderDimensions) {

        final boolean isEmpty = geometry.isEmpty();
        final GeometryType geometryType = GeometryType.valueOf(geometry.getClass());
        TWKBHeader header = forcePreserveHeaderDimensions ? params
            : setDimensions(geometry, params);
        header = header.setEmpty(isEmpty).setGeometryType(geometryType);

        if (optimizedEncoding) {
            if (isEmpty && header.hasExtendedPrecision()) {
                header = header.setHasZ(false).setHasM(false);
            }
            if ((isEmpty || geometryType == POINT) && header.hasBBOX()) {
                header = header.setHasBBOX(false);
            }
        } else {
            if (isEmpty && header.hasBBOX()) {
                header = header.setHasBBOX(false);
            }
        }
        return header;
    }

    private void writeGeometryBody(Geometry geom, TWKBOutputStream out, TWKBHeader header)
        throws IOException {
        if (header.isEmpty()) {
            return;
        }
        if (header.hasBBOX()) {
            writeBbox(geom, out, header);
        }
        final GeometryType geometryType = GeometryType.valueOf(geom.getClass());
        switch (geometryType) {
            case POINT:
                writePoint((Point) geom, out, header);
                return;
            case LINESTRING:
                writeLineString((LineString) geom, out, header, new long[header.getDimensions()]);
                return;
            case POLYGON:
                writePolygon((Polygon) geom, out, header, new long[header.getDimensions()]);
                return;
            case MULTIPOINT:
                writeMultiPoint((MultiPoint) geom, out, header);
                return;
            case MULTILINESTRING:
                writeMultiLineString((MultiLineString) geom, out, header);
                return;
            case MULTIPOLYGON:
                writeMultiPolygon((MultiPolygon) geom, out, header);
                return;
            case GEOMETRYCOLLECTION:
                writeGeometryCollection((GeometryCollection) geom, out, header);
                return;
            default:
                break;
        }
    }

    private void writePoint(Point geom, TWKBOutputStream out, TWKBHeader header)
        throws IOException {
        assert !geom.isEmpty();
        CoordinateSequence seq = geom.getCoordinateSequence();
        int dimensions = header.getDimensions();
        for (int d = 0; d < dimensions; d++) {
            writeOrdinate(seq.getOrdinate(0, d), 0L, header.getPrecision(d), out);
        }
    }

    private void writeCoordinateSequence(CoordinateSequence coordinateSequence,
        TWKBOutputStream out, TWKBHeader header, long[] prev) throws IOException {
        int size = coordinateSequence.size();
        out.writeUnsignedVarInt(size);
        writeCoordinateSequence(coordinateSequence, size, out, header, prev);
    }

    private void writeCoordinateSequence(CoordinateSequence coordinateSequence, int size,
        TWKBOutputStream out, TWKBHeader header, long[] prev) throws IOException {

        final int dimensions = header.getDimensions();
        for (int coordIndex = 0; coordIndex < size; coordIndex++) {
            for (int ordinateIndex = 0; ordinateIndex < dimensions; ordinateIndex++) {
                long previousValue = prev[ordinateIndex];
                int precision = header.getPrecision(ordinateIndex);
                double ordinate = coordinateSequence.getOrdinate(coordIndex, ordinateIndex);
                long preciseOrdinate = writeOrdinate(ordinate, previousValue, precision, out);
                prev[ordinateIndex] = preciseOrdinate;
            }
        }
    }

    private long writeOrdinate(double ordinate, long previousOrdinateValue, int precision,
        TWKBOutputStream out) throws IOException {
        long preciseOrdinate = makePrecise(ordinate, precision);
        long delta = preciseOrdinate - previousOrdinateValue;
        out.writeSignedVarLong(delta);
        return preciseOrdinate;
    }

    private long makePrecise(double value, int precision) {
        return Math.round(value * Math.pow(10, precision));
    }

    private void writeLineString(LineString geom, TWKBOutputStream out, TWKBHeader header,
        long[] prev) throws IOException {
        writeCoordinateSequence(geom.getCoordinateSequence(), out, header, prev);
    }

    private void writePolygon(Polygon geom, TWKBOutputStream out, TWKBHeader header,
        long[] prev) throws IOException {
        if (geom.isEmpty()) {
            out.writeUnsignedVarInt(0);
            return;
        }
        final int numInteriorRing = geom.getNumInteriorRing();
        final int nrings = 1 + numInteriorRing;
        out.writeUnsignedVarInt(nrings);
        writeLinearRing(geom.getExteriorRing(), out, header, prev);
        for (int r = 0; r < numInteriorRing; r++) {
            writeLinearRing(geom.getInteriorRingN(r), out, header, prev);
        }
    }

    private void writeLinearRing(LinearRing geom, TWKBOutputStream out, TWKBHeader header,
        long[] prev) throws IOException {
        if (geom.isEmpty()) {
            out.writeUnsignedVarInt(0);
            return;
        }
        CoordinateSequence seq = geom.getCoordinateSequence();
        int size = seq.size();
        if (optimizedEncoding && seq.size() > 2) {
            // With linear rings we can save one coordinate, they're automatically closed at parsing
            // time. But we can only do that if due to precision lost the two endpoints won't be
            // equal, otherwise the parser won't know it has to close the linear ring
            double x1 = seq.getOrdinate(0, 0);
            double y1 = seq.getOrdinate(0, 1);
            double x2 = seq.getOrdinate(size - 2, 0);
            double y2 = seq.getOrdinate(size - 2, 1);
            int precision = header.getPrecision(0);
            if (makePrecise(x1, precision) != makePrecise(x2, precision)
                || makePrecise(y1, precision) != makePrecise(y2, precision)) {
                --size;
            }
        }
        out.writeUnsignedVarInt(size);
        writeCoordinateSequence(seq, size, out, header, prev);
    }

    private void writeMultiPoint(MultiPoint geom, TWKBOutputStream out, TWKBHeader header)
        throws IOException {
        assert !geom.isEmpty();

        CoordinateSequence seq = geom.getFactory().getCoordinateSequenceFactory()
            .create(geom.getCoordinates());
        writeCoordinateSequence(seq, out, header, new long[header.getDimensions()]);
    }

    private void writeMultiLineString(MultiLineString geom, TWKBOutputStream out,
        TWKBHeader header) throws IOException {
        final int size = writeNumGeometries(geom, out);
        long[] prev = new long[header.getDimensions()];
        for (int i = 0; i < size; i++) {
            writeLineString((LineString) geom.getGeometryN(i), out, header, prev);
        }
    }

    private void writeMultiPolygon(MultiPolygon geom, TWKBOutputStream out,
        TWKBHeader header) throws IOException {
        final int size = writeNumGeometries(geom, out);
        long[] prev = new long[header.getDimensions()];
        for (int i = 0; i < size; i++) {
            writePolygon((Polygon) geom.getGeometryN(i), out, header, prev);
        }
    }

    private void writeGeometryCollection(GeometryCollection geom, TWKBOutputStream out,
        TWKBHeader header) throws IOException {
        final int size = writeNumGeometries(geom, out);
        for (int i = 0; i < size; i++) {
            Geometry geometryN = geom.getGeometryN(i);
            boolean forcePreserveDimensions = geometryN.isEmpty();
            write(geometryN, out, header, forcePreserveDimensions);
        }
    }

    private int writeNumGeometries(GeometryCollection geom, TWKBOutputStream out)
        throws IOException {
        int size = geom.getNumGeometries();
        out.writeUnsignedVarInt(size);
        return size;
    }

    private void writeBbox(Geometry geom, TWKBOutputStream out, TWKBHeader header)
        throws IOException {
        final int dimensions = header.getDimensions();
        final double[] boundsCoordinates = computeEnvelope(geom, dimensions);

        for (int d = 0; d < dimensions; d++) {
            final int precision = header.getPrecision(d);
            double min = boundsCoordinates[2 * d];
            double max = boundsCoordinates[2 * d + 1];
            long preciseMin = writeOrdinate(min, 0, precision, out);
            writeOrdinate(max, preciseMin, precision, out);
        }
    }

    private static double[] computeEnvelope(Geometry geom, int dimensions) {
        BoundsExtractor extractor = new BoundsExtractor(dimensions);
        geom.apply(extractor);
        return extractor.ordinates;
    }

    private static TWKBHeader setDimensions(Geometry g, TWKBHeader header) {
        if (g.isEmpty()) {
            return header.setHasZ(false).setHasM(false);
        }
        if (g instanceof Point) {
            return setDimensions(((Point) g).getCoordinateSequence(), header);
        }
        if (g instanceof LineString) {
            return setDimensions(((LineString) g).getCoordinateSequence(), header);
        }
        if (g instanceof Polygon) {
            return setDimensions(((Polygon) g).getExteriorRing().getCoordinateSequence(), header);
        }
        return setDimensions(g.getGeometryN(0), header);
    }

    private static TWKBHeader setDimensions(CoordinateSequence seq, TWKBHeader header) {
        boolean hasZ = seq.hasZ();
        boolean hasM = seq.hasM();
        return header.setHasZ(hasZ).setHasM(hasM);
    }

    private static class BoundsExtractor implements CoordinateSequenceFilter {

        private final boolean done = false;

        private final boolean geometryChanged = false;

        private final int dimensions;

        double[] ordinates = new double[] { //
            Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, // note, Double.MIN_VALUE is positive
            Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, //
            Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, //
            Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY//
        };

        BoundsExtractor(int dimensions) {
            this.dimensions = dimensions;
        }

        public @Override void filter(final CoordinateSequence seq, final int coordIndex) {
            for (int ordinateIndex = 0; ordinateIndex < dimensions; ordinateIndex++) {
                final double ordinate = seq.getOrdinate(coordIndex, ordinateIndex);
                final int minIndex = 2 * ordinateIndex;
                final int maxIndex = minIndex + 1;
                double minValue = ordinates[minIndex];
                double maxValue = ordinates[maxIndex];
                minValue = Math.min(minValue, ordinate);
                maxValue = ordinate > maxValue ? ordinate : maxValue;// Math.max(maxValue, ordinate);
                ordinates[minIndex] = minValue;
                ordinates[maxIndex] = maxValue;
            }
        }

        @Override
        public boolean isDone() {
            return this.done;
        }

        @Override
        public boolean isGeometryChanged() {
            return this.geometryChanged;
        }
    }
}