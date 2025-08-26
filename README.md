Rotate some ocean bottom seismo data with me!
Still a work in progress getting the flow down but heres a short list.

1. mass download your data using what ever you like

2. verify that miniseed files exist 

3. in scripts/miniseed_sac.py replace event_id with the name of the event/directory. Run scripts/miniseed_sac.py. you should now have sac files in the sac directory

4. replace event name in scripts/update_xml.py. run this script

5. update event name in 2_12Z_RTZ.py and run 

rotated .sac files are now sac_rotated. you can move these to the sac directory if you'd like or keep em there. up to you!

working on the work flow right now to make it automated for lots of events but this is a good time to create a repo for it. 
