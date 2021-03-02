# Orquestra Resource To Run Madness with Tequila
- uses the docker container:  [kottmanj/tequila-madness](https://dockerhub.com/kottmanj/tequila-madness)
- example workflow is provided in [example_workflow](example_workflow)

# Extending/Improving this
- make a fork of this repo
- adapt the github address in the example_workflow accordingly (so that it is now your fork)
- later make a pull request

# Using this

Execute like this
```bash
qe submit workflow my_workflow.yaml
```

Get results
```bash
qe get workflowresults WORKFLOW-ID
```

Debug like this  
```bash
# get information about STEP-ID and which step crashed
qe get workflow WORKFLOW-ID
# get error logs
qe get logs WORKFLOW-ID -s STEP-ID
```

See Orquestra [docs](http://docs.orquestra.io/) for more

# TODO
- need to have tequila molecules (at least the madness ones) JSON serializable (so that we can dump and load)
- make some example workflows (run madness and then do several tq tasks in parallel with the obtained molecule)
- figure out how to make own `language` in orquestra, so that we can use the tequila installed in the docker container instead of re-installing (easier version control)

# Dockerfile
[Here](Docker/Dockerfile) you can find the Dockerfile that produces the image hosten on [kottmanj/tequila-madness](https://dockerhub.com/kottmanj/tequila-madness).  
You don't need this file, but you can modify it and host your own image. Works like this:
1. Modify the Dockerfile
2. Make a docker image (takes a while)
```bash
sudo docker image -t username/imagename:version Dockerfile
```
3. Push the image to dockerhub
```bash
sudo docker push username/imagename:version
```

The Dockerfile currently uses ubuntu. Starting from another base might be beneficial for lighter containes.  
For now we are happy with this.  
