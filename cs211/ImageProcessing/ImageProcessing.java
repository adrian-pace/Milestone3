package cs211.ImageProcessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import processing.core.PApplet;
import processing.core.PImage;
import processing.core.PVector;

public class ImageProcessing extends PApplet {
	PImage img;
	HScrollbar thresholdBar1;
	HScrollbar thresholdBar2;
	HScrollbar thresholdBar3;
	HScrollbar thresholdBar4;

	public void setup() {
		size(800, 680);
		img = loadImage("board1.jpg");
		thresholdBar1 = new HScrollbar(this, 0, 600, 800, 20);
		thresholdBar2 = new HScrollbar(this, 0, 620, 800, 20);
		thresholdBar3 = new HScrollbar(this, 0, 640, 800, 20);
		thresholdBar4 = new HScrollbar(this, 0, 660, 800, 20);

		// noLoop(); // no interactive behaviour: draw() will be called only
		// once.
	}

	public void draw() {
		background(color(0, 0, 0));
		// ...
		int lowThreshold = (int) (thresholdBar1.getPos() * 255);
		int highThreshold = (int) (thresholdBar2.getPos() * 255);

		/* if only one convolute, put hue to 116,136. */
		PImage result = hueImage(img, 115, 135);
		result = convolute(result);
		// result=thresholdBinary(result, 30, false);

		result = sobel(result, 0.17);
		// PImage houghImage = hough(result);
		image(result, 0, 0);
		// image(result, 0, 0);
		List<PVector> lines = hough(result, 6);
		QuadGraph quadGraph = new QuadGraph();
		quadGraph.build(lines, width, height);
		quadGraph.findCycles();

		for (int[] quad : quadGraph.cycles) {
			PVector l1 = lines.get(quad[0]);
			PVector l2 = lines.get(quad[1]);
			PVector l3 = lines.get(quad[2]);
			PVector l4 = lines.get(quad[3]);
			// (intersection() is a simplified version of the
			// intersections() method you wrote last week, that simply
			// return the coordinates of the intersection between 2 lines)

			PVector c12 = intersection(l1, l2);
			PVector c23 = intersection(l2, l3);
			PVector c34 = intersection(l3, l4);
			PVector c41 = intersection(l4, l1);
			if (QuadGraph.isConvex(c12, c23, c34, c41)
					&& QuadGraph.nonFlatQuad(c12, c23, c34, c41)
					&& QuadGraph.validArea(c12, c23, c34, c41, 400000,
							20000)) {
				// Choose a random, semi-transparent colour
				Random random = new Random();
				fill(color(min(255, random.nextInt(300)),
						min(255, random.nextInt(300)),
						min(255, random.nextInt(300)), 50));
				quad(c12.x, c12.y, c23.x, c23.y, c34.x, c34.y, c41.x, c41.y);
			}
		}

		List<PVector> intersections = getIntersections(lines);

		// thresholdBar1.display();
		// thresholdBar1.update();
		// thresholdBar2.display();
		// thresholdBar2.update();
		// thresholdBar3.display();
		// thresholdBar3.update();
		// thresholdBar4.display();
		// thresholdBar4.update();
		// System.out.println(thresholdBar1.getPos() * 255 + " "
		// + thresholdBar2.getPos() * 255);
		// System.out.println(thresholdBar4.getPos());
	}

	private PVector intersection(PVector line1, PVector line2) {
		double d = Math.cos(line2.y) * Math.sin(line1.y) - Math.cos(line1.y)
				* Math.sin(line2.y);
		int x = (int) Math.round((line2.x * Math.sin(line1.y) - line1.x
				* Math.sin(line2.y))
				/ d);
		int y = (int) Math.round((-line2.x * Math.cos(line1.y) + line1.x
				* Math.cos(line2.y))
				/ d);
		return new PVector(x, y);
	}

	public ArrayList<PVector> hough(PImage edgeImg, int nLines) {
		int threshold = 200;
		float discretizationStepsPhi = 0.06f;
		float discretizationStepsR = 2.5f;
		// dimensions of the accumulator
		int phiDim = (int) (Math.PI / discretizationStepsPhi);
		int rDim = (int) (((edgeImg.width + edgeImg.height) * 2 + 1) / discretizationStepsR);
		// our accumulator (with a 1 pix margin around)
		int[] accumulator = new int[(phiDim + 2) * (rDim + 2)];
		// Fill the accumulator: on edge points (ie, white pixels of the edge
		// image), store all possible (r, phi) pairs describing lines going
		// through the point.
		for (int y = 0; y < edgeImg.height; y++) {
			for (int x = 0; x < edgeImg.width; x++) {
				// Are we on an edge?
				if (brightness(edgeImg.pixels[y * edgeImg.width + x]) != 0) {
					// ...determine here all the lines (r, phi) passing through
					// pixel (x,y), convert (r,phi) to coordinates in the
					// accumulator, and increment accordingly the accumulator.
					for (int i = 0; i < phiDim; i++) {
						double phi = i * discretizationStepsPhi;
						double r = x * Math.cos(phi) + y * Math.sin(phi);
						int accPhi = i;
						int accR = (int) (r / discretizationStepsR + (rDim - 1) / 2);
						accumulator[(accPhi + 1) * (rDim + 2) + accR + 1]++;
					}

				}
			}
		}
		PImage houghImg = createImage(rDim + 2, phiDim + 2, ALPHA);
		for (int i = 0; i < accumulator.length; i++) {
			houghImg.pixels[i] = color(min(255, accumulator[i]));
		}
		houghImg.updatePixels();
		// return houghImg;

		ArrayList<Integer> bestCandidates = new ArrayList<Integer>();
		// size of the region we search for a local maximum
		int neighbourhood = 10;
		// only search around lines with more that this amount of votes
		// (to be adapted to your image)
		int minVotes = 200;
		for (int accR = 0; accR < rDim; accR++) {
			for (int accPhi = 0; accPhi < phiDim; accPhi++) {
				// compute current index in the accumulator
				int idx = (accPhi + 1) * (rDim + 2) + accR + 1;
				if (accumulator[idx] > minVotes) {
					boolean bestCandidate = true;
					// iterate over the neighbourhood
					for (int dPhi = -neighbourhood / 2; dPhi < neighbourhood / 2 + 1; dPhi++) {
						// check we are not outside the image
						if (accPhi + dPhi < 0 || accPhi + dPhi >= phiDim)
							continue;
						for (int dR = -neighbourhood / 2; dR < neighbourhood / 2 + 1; dR++) {
							// check we are not outside the image
							if (accR + dR < 0 || accR + dR >= rDim)
								continue;
							int neighbourIdx = (accPhi + dPhi + 1) * (rDim + 2)
									+ accR + dR + 1;
							if (accumulator[idx] < accumulator[neighbourIdx]) {
								// the current idx is not a local maximum!
								bestCandidate = false;
								break;
							}
						}
						if (!bestCandidate)
							break;
					}
					if (bestCandidate) {
						// the current idx *is* a local maximum
						bestCandidates.add(idx);
					}
				}
			}
		}

		Collections.sort(bestCandidates, new HoughComparator(accumulator));
		ArrayList<PVector> returnArray = new ArrayList<PVector>();
		for (int i = 0; i < nLines && i < bestCandidates.size(); i++) {

			// first, compute back the (r, phi) polar coordinates:
			int accPhi = (int) (bestCandidates.get(i) / (rDim + 2)) - 1;
			int accR = bestCandidates.get(i) - (accPhi + 1) * (rDim + 2) - 1;
			float r = (accR - (rDim - 1) * 0.5f) * discretizationStepsR;
			float phi = accPhi * discretizationStepsPhi;
			returnArray.add(new PVector(r, phi));

			// Cartesian equation of a line: y = ax + b
			// in polar, y = (-cos(phi)/sin(phi))x + (r/sin(phi))
			// => y = 0 : x = r / cos(phi)
			// => x = 0 : y = r / sin(phi)
			// compute the intersection of this line with the 4 borders of
			// the image
			int x0 = 0;
			int y0 = (int) (r / sin(phi));
			int x1 = (int) (r / cos(phi));
			int y1 = 0;
			int x2 = edgeImg.width;
			int y2 = (int) (-cos(phi) / sin(phi) * x2 + r / sin(phi));
			int y3 = edgeImg.width;
			int x3 = (int) (-(y3 - r / sin(phi)) * (sin(phi) / cos(phi)));
			// Finally, plot the lines
			stroke(204, 102, 0);
			if (y0 > 0) {
				if (x1 > 0)
					line(x0, y0, x1, y1);
				else if (y2 > 0)
					line(x0, y0, x2, y2);
				else
					line(x0, y0, x3, y3);
			} else {
				if (x1 > 0) {
					if (y2 > 0)
						line(x1, y1, x2, y2);
					else
						line(x1, y1, x3, y3);
				} else
					line(x2, y2, x3, y3);
			}

		}
		return returnArray;
	}

	ArrayList<PVector> getIntersections(List<PVector> lines) {
		ArrayList<PVector> intersections = new ArrayList<PVector>();
		for (int i = 0; i < lines.size() - 1; i++) {
			PVector line1 = lines.get(i);
			for (int j = i + 1; j < lines.size(); j++) {
				PVector line2 = lines.get(j);
				// compute the intersection and add it to 'intersections'
				double d = Math.cos(line2.y) * Math.sin(line1.y)
						- Math.cos(line1.y) * Math.sin(line2.y);
				int x = (int) Math.round((line2.x * Math.sin(line1.y) - line1.x
						* Math.sin(line2.y))
						/ d);
				int y = (int) Math
						.round((-line2.x * Math.cos(line1.y) + line1.x
								* Math.cos(line2.y))
								/ d);
				intersections.add(new PVector(x, y));
				// draw the intersection
				fill(255, 128, 0);
				ellipse(x, y, 10, 10);
			}
		}
		return intersections;
	}

	PImage sobel(PImage img, double d) {
		int[][] hKernel = { { 0, 1, 0 }, { 0, 0, 0 }, { 0, -1, 0 } };
		int[][] vKernel = { { 0, 0, 0 }, { 1, 0, -1 }, { 0, 0, 0 } };
		PImage result = createImage(img.width, img.height, ALPHA);
		result.loadPixels();
		img.loadPixels();
		int dim = img.width * img.height;
		float sumh = 0;
		float sumv = 0;
		float[] buffer = new float[dim];
		float max = 0;
		for (int y = 1; y < img.height - 1; y++) {
			for (int x = 1; x < img.width - 1; x++) {
				sumh = 0;
				sumv = 0;
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						sumh += brightness(img.pixels[(y * img.width + x + j
								* img.width + k + dim)
								% dim])
								* hKernel[j + 1][k + 1];
						sumv += brightness(img.pixels[(y * img.width + x + j
								* img.width + k + dim)
								% dim])
								* vKernel[j + 1][k + 1];
						buffer[y * img.width + x] = sqrt(sumh * sumh + sumv
								* sumv);
						if (buffer[y * img.width + x] > max)
							max = buffer[y * img.width + x];
					}

				}
				// result.pixels[i]=;
			}

		}
		for (int y = 2; y < img.height - 2; y++) { // Skip top and bottom edges
			for (int x = 2; x < img.width - 2; x++) { // Skip left and right
				if (buffer[y * img.width + x] > (int) (max * d)) {
					result.pixels[y * img.width + x] = color(255);
				} else {
					result.pixels[y * img.width + x] = color(0);
				}
			}
		}
		result.updatePixels();
		return result;
	}

	PImage convolute(PImage img) {
		int[][] matrix = { { 9, 12, 9 }, { 12, 15, 12 }, { 9, 12, 9 } };
		PImage result = createImage(img.width, img.height, ALPHA);
		result.loadPixels();
		img.loadPixels();
		int dim = img.width * img.height;
		int weight = 0;
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				weight += matrix[i][j];
			}
		}
		int sum = 0;
		for (int y = 1; y < img.height - 1; y++) {
			for (int x = 1; x < img.width - 1; x++) {
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						sum += brightness(img.pixels[(y * img.width + x + j
								* img.width + k + dim)
								% dim])
								* matrix[j + 1][k + 1];
					}
				}
				result.pixels[x + y * img.width] = color(sum / weight);
				sum = 0;
			}
		}
		result.updatePixels();
		return result;
	}

	PImage hueImage(PImage img, int lowThreshold, int highThreshold) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			if (hue(img.pixels[i]) <= highThreshold
					&& hue(img.pixels[i]) >= lowThreshold)
				result.pixels[i] = img.pixels[i];
			else
				result.pixels[i] = color(0);
		}
		result.updatePixels();
		return result;
	}

	PImage hueImage(PImage img) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			result.pixels[i] = color((int) hue(img.pixels[i]));
		}
		result.updatePixels();
		return result;
	}

	PImage thresholdBinary(PImage img, int threshold, boolean inverted) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			if (inverted) {
				if (brightness(img.pixels[i]) < threshold) {
					result.pixels[i] = img.pixels[i];
				} else {
					result.pixels[i] = color(0);
				}
			} else {
				if (brightness(img.pixels[i]) > threshold) {
					result.pixels[i] = img.pixels[i];
				} else {
					result.pixels[i] = color(0);
				}

			}

		}
		result.updatePixels();
		return result;

	}
}
