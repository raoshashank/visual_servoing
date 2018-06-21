import rospy
from std_msgs.msg import Int64
import cv2.cv2
import cv2 as cv
import numpy as np
from object_tracking.msg import coordinates
def main():
    global pub,coo,cap,cX,cY
    coo.x=0
    coo.y=0
    cap=cv2.VideoCapture(0)
    while not rospy.is_shutdown():
       while (cap.isOpened()):
        ret, frame = cap.read()
        # cv2.imshow('frame',frame)

        hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        lower_bound = np.array([25, 50, 50])
        upper_bound = np.array([50, 255, 255])
        mask = cv2.inRange(hsv, lower_bound, upper_bound)
        res = cv2.bitwise_and(frame, frame, mask=mask)
        kernel_e = np.ones((3, 3), dtype=np.uint8)
        res2 = cv2.erode(res, kernel=kernel_e, iterations=2)
        res3 = cv2.dilate(res2, kernel=kernel_e, iterations=2)
        # blur=cv.GaussianBlur(res3,(5,5),0)
        blur_gray = cv2.cvtColor(res3, cv2.COLOR_BGR2GRAY)
        _, cnts, _ = cv2.findContours(blur_gray, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        x=0
        y=0
        radius=0
        for i in cnts:
            if len(i) > 50:
                #M = cv.moments(i)
                #cX = int(M["m10"] / M["m00"])
                #cY = int(M["m01"] / M["m00"])
                (x,y),radius = cv2.minEnclosingCircle(i)
                #cv2.drawContours(res3, [i], -1, (0, 255, 0), 2)
                #cv2.circle(res3, (cX, cY), 7, (255, 255, 255), -1)
                #cv2.putText(res3, "center", (cX - 20, cY - 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 2)
  
        coo.x=x
        coo.y=y
        coo.radius=radius
        pub.publish(coo) 
              
             
                # cv2.drawContours(res3, [i], -1, (0, 255, 0), 2)
                # cv2.circle(res3, (cX, cY), 7, (255, 255, 255), -1)
                # cv2.putText(res3, "center", (cX - 20, cY - 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 2)
        # cv.imshow('at',at)
        # cv2.imshow('res',res)
        # cv2.imshow('res2',res2)
        cv2.imshow('res3',res3)
        # cv2.imshow('blur', blur_gray)
        # cv2.imshow('res3', res3)
        # # cv.imshow('gray', gray)
        # cv2.imshow('mask',mask)
        # cv2.imshow('hsv',hsv)
    
        if cv2.waitKey(1) & 0xFF == ord('q'):
            cap.release()
            cv2.destroyAllWindows()

if __name__=="__main__":
    global pub,coo,cap,cX,cY
    cX=0
    cY=0
    rospy.init_node('tracker',anonymous=False) 
    pub=rospy.Publisher('/this_topic',coordinates,queue_size=1)
    coo=coordinates()
    main()
    #rospy.spin()
