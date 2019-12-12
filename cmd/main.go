/*
 *  main.go
 *  cmd
 *
 *  Created by Haibao Tang on 12/12/19
 *  Copyright © 2019 Haibao Tang. All rights reserved.
 */

package main

import (
	"github.com/op/go-logging"
	"github.com/tanghaibao/allhic"
)

// main is the entrypoint for the entire program, routes to commands
func main() {
	logging.SetBackend(allhic.BackendFormatter)
	allhic.Execute()
}
