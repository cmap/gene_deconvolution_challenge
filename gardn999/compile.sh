mkdir -p build
javac src/*.java -d build
cd build

jar cfe ../gardn999_CMap_DPeak.jar Main *.class

cd ..
echo "Jar file successfully created."
