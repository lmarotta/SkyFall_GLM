%training data
load sensors_data

boxplot(acce);
title('Training Accelerometer Data')

boxplot(gyro);
title('Training Gyroscope Data')

boxplot(baro);
title('Training Barometer Data')

%phone data
load accTest
acce = [X Y Z]
load gyrTest
gyro = [X Y Z]
load barTest
baro = [Altitude Pressure]

boxplot(acce);
title('Phone Accelerometer Data')

boxplot(gyro);
title('Phone Gyroscope Data')

boxplot(baro);
title('Phone Barometer Data')