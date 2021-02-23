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
