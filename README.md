### Robot Motion Planning in dimensions greater than 3D.

#### Implementation
The implementation is in Python and from scratch.

All obstacles and the boundary are modeled using ellipsoids as modeling of real-world obstacles as ellipsoids is a solved problem.

#### Examples

#### Input
The input files are named as such. The general formulation is:
First line : the first integer is the number of obstacles and the second integer is the number of dimensions

#### Motivation
The Silhouette method is complete. Hence, finding a path using the same and then optimizing it helps ensure completeness of the planning method. It improves on methods that use Heuristics to find a path, but lack a guarantee of completeness.

#### Installation
Clone this repo, it uses python 2.7 and its associate libraries.

#### Performance
The Silhouette method is complete. Hence if a path exists from the start to the goal point, it will find it.

#### Usability
This project contains a python based GUI that is also portable, such that this software can run on any system, even one that do not have a python installation.

#### The MIT License
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
