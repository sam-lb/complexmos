



class Node {

    static nodeID = 0;

    constructor(value) {
        this.value = value;
        this.id = Node.nodeID;
        Node.nodeID++;
        this.next = null;   // horizontal child
        this.child = null;  // vertical child
    }

    setNext(next) {
        this.next = next;
    }

    hasNext() {
        return this.next !== null;
    }

    toString() {
        return `Node(${this.value})`;
    }

}


class LinkedList {

    constructor(valuesArray=null) {
        this.head = null;
        this.size = 0;

        if (valuesArray !== null) {
            for (let value of valuesArray) {
                this.push(value);
            }
        }
    }

    push(value) {
        const newNode = new Node(value);
        
        if (this.head === null) {
            this.head = newNode;
        } else {
            let node = this.head;
            while (node.hasNext()) {
                node = node.next;
            }
            node.setNext(newNode);
        }
        this.size++;
    }

    remove(value) {
        if (this.head === null) {
            return -1;
        } else {
            if (this.head.value === value) {
                this.head = this.head.next;
                this.size--;
                return 1;
            } else {
                let node = this.head;
                while (node.hasNext()) {
                    if (node.next.value === value) {
                        node.next = node.next.next;
                        this.size--;
                        return 1;
                    }
                    node = node.next;
                }
                return -1;
            }
        }
    }

    toString() {        
        if (this.head === null) {
            return "LinkedList([])";
        } else {
            let result = "[";
            let node = this.head;
            while (node !== null) {
                result += node.toString() + ", ";
                node = node.next;
            }
            result = result.slice(0, -2) + "]";
            return result;
        }
    }

}