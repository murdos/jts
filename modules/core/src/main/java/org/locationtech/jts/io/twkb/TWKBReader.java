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

import java.io.ByteArrayInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.impl.PackedCoordinateSequenceFactory;
import org.locationtech.jts.io.ParseException;

public class TWKBReader {

    private static final GeometryFactory DEFAULT_FACTORY = new GeometryFactory(
        PackedCoordinateSequenceFactory.DOUBLE_FACTORY);

    private GeometryFactory geometryFactory;

    public TWKBReader() {
        this(DEFAULT_FACTORY);
    }

    public TWKBReader(GeometryFactory geometryFactory) {
        this.geometryFactory = geometryFactory;
    }

    public TWKBReader setGeometryFactory(GeometryFactory geometryFactory) {
        this.geometryFactory = geometryFactory;
        return this;
    }

    public Geometry read(byte[] bytes) throws ParseException {
        return read(new ByteArrayInputStream(bytes));
    }

    public Geometry read(InputStream in) throws ParseException {
        return read((DataInput) new DataInputStream(in));
    }

    public Geometry read(DataInput in) throws ParseException {
        try {
            return TWKBIO.read(geometryFactory, in);
        } catch (IOException ex) {
            throw new ParseException("Unexpected IOException caught: " + ex.getMessage());
        }
    }
}
